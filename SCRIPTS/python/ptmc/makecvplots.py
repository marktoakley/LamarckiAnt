import matplotlib as mpl
mpl.use("Agg") # so we can use it without an X server

import itertools
import os
import getopt, sys
import copy
import numpy as np
import matplotlib.pyplot as plt

import histogram_reweighting.histogram_reweighting1d as WHAM

import ptmc_utils as utils


def makeHistogram(Vall, index, indexsub=None, neglect=0.01, strip=True, verbose=False):
    """
    generate a histogram from a list of VisitsHis objects
    
    Parameters
    ----------
    Vall : list of lists of VisitsHis objects
        Vall[r][i] is the VisitsHis object of replica r at time index i
    index : int
        the time index for which to create the histogram
    indexsub : int, optional
        the data at time index indexsub will be subtracted so only the data from
        time interval indexsub to index will be included
    neglect : float, optional
        neglect bins which have fewer visits than neglect*max_visits, where max_visits
        is the maximum visits for that replica 
    strip : bool, optional
        remove all bins from the front and back that are empty for all replicas
    """
    vis = Vall[0][index]
    vis.read()
    binenergy = np.copy(vis.energies_all)
    nebins = len(binenergy)
    nreps = len(Vall)
    visits = np.zeros( [nreps, nebins] )
    if verbose:
        print "\n\n"
        print "#######################################################################"
        if indexsub != None:
            vis2 = Vall[0][indexsub]
            print "histogram", Vall[0][index].fname, "minus histogram", vis2.fname
        else:
            print "histogram", Vall[0][index].fname
        print "#######################################################################"


    for k in range(nreps):
        vis = Vall[k][index]
        vis.read()


        if indexsub != None:
            vis2 = Vall[k][indexsub]
            vis2.read()
            e, counts = vis.getHisDiffAll( vis2 )
            #print counts
            #enew, countsnew = vis.getHisAll()
            #enew2, countsnew2 = vis2.getHisAll()
            #np.savetxt("c1", countsnew)
            #np.savetxt("c2", countsnew2)
            #np.savetxt("c3", countsnew-countsnew2)
            #print countsnew - countsnew2
            #exit()
        else:
            counts = vis.counts_all
        visits[k,:] = counts
        
    return prepare_visits(binenergy, visits, neglect=neglect, strip=strip, verbose=verbose)

def prepare_visits(binenergy, visits, neglect=None, strip=True, verbose=False):
    """clean up a visits collection
    
    Parameters
    ----------
    binenergy : ndarray
        the (lower bound) energy of the bins
    visits : 2d array
        visits[r,b] is the number of visits in replica r to energy bin b
    neglect : float
        if there are too few visits in a bin, set it to zero.  For a given replica the 
        cutoff is neglect times the maximum number of visits to a bin
    strip : bool 
        remove all bins from the front and back that are empty for all replicas

    """
    #if there are too few visits in a bin, set it to zero
    if neglect is not None:
        nreps = visits.shape[0]
        for k in range(nreps):
            visitsmax=visits[k,:].max()
            neglectmax = float(visitsmax) * neglect
            neglectmax = max(neglectmax, 1)
            neglectmax = min(neglectmax, 10000)
            ind=np.where(visits[k,:] <= neglectmax)
            if verbose:
                print k, "total visits ", visits[k,:].sum(), "max visits ", visitsmax, "max visits to be discarded ", 
                print visits[k,ind[0]].max(), visits[k,ind[0]].max()/ visitsmax #, "bins discarded", len(ind[0]), "out of", nnonzero
            visits[k,ind[0]] = 0


    if strip:
        #remove all bins from the front and back that are empty for all replicas
        nonzero = np.where( visits.sum(0) != 0 )[0]
        ibegin = nonzero[0]
        iend = nonzero[-1]+1  #one past the last index
        binenergy = binenergy[ibegin : iend]
        visits = visits[:,ibegin : iend]
    
    #np.savetxt("c4", visits)
    return binenergy, visits

def getTlist(fname = "temperatures"):
    """try to find the temperatures of the replicas by reading file fname"""
    Tlist = []
    with open(fname, "r") as fin:
        for line in fin:
            Tlist.append( float( line.split()[0] ) )
    return Tlist

def checkWHAM(binenergy, visits, cvdata, Tlist):
    for i in range(len(Tlist)):
        m, std = utils.wstd(binenergy, visits[i,:])
        print "wham_check", Tlist[i], m, std

def calcCvNoWHAM(binenergy, visits, Tlist, NDOF):
    """compute the heat capacity using the direct calculation, i.e. not WHAM
    
    This method computes <E^2> - <E>^2 directly from the visits histogram.
    The heat capacity is returned at the temperatures of the replicas.  
    
    Parameters
    ----------
    binenergy : ndarray
        the (lower bound) energy of the bins
    visits : 2d array
        visits[r,b] is the number of visits in replica r to energy bin b
    Tlist : array
        the temperatures of the replicas
    ndof : int
        the number of degrees of freedom in the system
    """
    nreps, nbins = np.shape(visits)
    cvdata = np.zeros([nreps, 6])
    for i in range(nreps):
        cvdata[i,0] = Tlist[i]
        m, std = utils.wstd(binenergy, visits[i,:])
        cvdata[i,4] = m
        cvdata[i,3] = std
        #print std, Tlist[i], NDOF, std**2, Tlist[i]**2
        cv = float(NDOF)/2. + std**2 / Tlist[i]**2
        cvdata[i,5] = cv
    return cvdata

def make_histogram_jackknife(Vcollection, isub, neglect=0.01):
    """make a histogram for the jackknife calculation
    
    
    Vcollection : 
        the collection of VisitsHis objects for each replica
    isub : int
        the index to subtract for the jackknife calculation

    Returns
    -------
    binenergy : float
        energy of the bins
    visits : 2d array
        the 2d array holding the visits histogram for each replica.  The data 
        collected during the time period between isub-1 and isub will
        be removed
    
    """
    nreps = Vcollection.get_nreplicas()
    
    visits = []
    for r in range(nreps):
        vf = Vcollection.vseries_list[r][-1]
        vi = Vcollection.vseries_list[r][0]
        v2 = Vcollection.vseries_list[r][isub]
        v1 = Vcollection.vseries_list[r][isub-1]
        
        vf.read()
        vi.read()
        v1.read()
        v2.read()
        
        counts = vf.counts_all - v2.counts_all + v1.counts_all - vi.counts_all
        visits.append(counts)
        
    binenergy = Vcollection.vseries_list[0][0].energies_all
    
    visits = np.array(visits)
    assert visits.shape == (Vcollection.get_nreplicas(), len(binenergy))
    
    return prepare_visits(binenergy, visits, neglect=neglect)
    

def estimate_error_jackknife(Vcollection, Tlist, NDOF=0, neglect=None):
    """estimate the error in Cv from the jackknife method
    
    The we have N data points X_1 ... X_N.  The error can be estimated
    by resampling the data.  We create N-1 sets of data X_[i] by taking
    data point X_i.  X_[i] = {X_k : k ! i}
    
    For each of these reduced data sets we compute the heat capacity.  Let sigma_cv
    be the square root of the variance of the heat capacity, then 
    sqrt(N-1) * sigma_cv is an unbiased estimator of the uncertainty in the
    heat capacity.
    
    In this case a data point X_i is the data (visits) collected between time 
    step i and time step i-1.  
    
    http://physics.ucsc.edu/~peter/jackboot.pdf
    """
    nvisits = Vcollection.get_nvisits()
    if nvisits <= 2:
        raise Exception("can't compute the error unless"
                        " the number of visits files is > 2: %d" % nvisits)
    
    heat_capacities = []
    weights = []
    
    total_time = Vcollection.get_time(-1) - Vcollection.get_time(0)
    
    for isub in range(1,nvisits):
        binenergy, visits = make_histogram_jackknife(Vcollection, isub, neglect=neglect)
        cvdata = calcCvNoWHAM(binenergy, visits, Tlist, NDOF)

        heat_capacities.append(cvdata[:,5])
        
        
        time2 = Vcollection.get_time(isub)
        time1 = Vcollection.get_time(isub-1)
        weights.append( total_time - (time2 - time1))
    
    ndata_sets = len(weights)
    
    heat_capacities = np.array(heat_capacities)
    weights = np.array(weights)
    # heat_capacities[i,t] is the heat capacity from dataset i at temperature t 
    assert heat_capacities.shape[1] == len(Tlist)
    
    # compute the mean and standard deviation of the heat capacites from each block
    cvmean = np.average(heat_capacities, weights=weights, axis=0)
    cv2mean = np.average(heat_capacities**2, weights=weights, axis=0)
#    cvvar = np.average( (heat_capacities - cvmean[np.newaxis,:])**2, weights=weights, axis=0 )
    cvvar = cv2mean - cvmean**2
#    cvvar *= float(ndata_sets) / (ndata_sets - 1)
#    cvstd = np.std(heat_capacities, axis=0)

    # we need to multiply by sqrt(ndata_sets) to get the unbiased estimator of the error
    cverr = np.sqrt(float(ndata_sets) * cvvar)
    
    
    # 
    # the above estimate of cvmean is biased, we can do better by using all the data at once
    #
    cvmean_jack = cvmean
    imax = Vcollection.get_nvisits() - 1
    isub = 0
    binenergy, visits = makeHistogram(Vcollection.vseries_list, imax, indexsub=isub, neglect=neglect)
    cvdata = calcCvNoWHAM(binenergy, visits, Tlist, NDOF)
    cvmean = cvdata[:,5]
    
    print "max difference between standard and jackkife Cv estimate", np.max(cvmean-cvmean_jack)

    
    return cvmean, cverr



#def estimate_error_bootstrap(Vcollection, Tlist, NDOF=0):
#    """estimate the error in Cv from the bootstrap method
#    
#    this is not the correct way to do the bootstrapping method.  I'm supposed
#    to do a random resampling of the data.
#    http://physics.ucsc.edu/~peter/jackboot.pdf
#    """
#    raise Exception("boostrap method is not implemented correctly")
#    nvisits = Vcollection.get_nvisits()
#    print "number of samples", nvisits
#    ndata_sets = nvisits - 1
#    
#    heat_capacities = []
#    weights = []
#    
#    for ind in range(1,nvisits):
#        isub = ind - 1
#
#        binenergy, visits = makeHistogram(Vcollection.vseries_list, ind, isub)
#        cvdata = calcCvNoWHAM(binenergy, visits, Tlist, NDOF)
#        heat_capacities.append(cvdata[:,5])
#        
#        time1 = Vcollection.get_time(ind)
#        time0 = Vcollection.get_time(isub)
#        weights.append( time1 - time0)
#    
#    heat_capacities = np.array(heat_capacities)
#    weights = np.array(weights)
#    # heat_capacities[i,t] is the heat capacity from dataset i at temperature t 
#    assert heat_capacities.shape[1] == len(Tlist)
#    
#    # compute the mean and standard deviation of the heat capacites from each block
#    cvmean = np.average(heat_capacities, weights=weights, axis=0)
#    cvvar = np.average( (heat_capacities - cvmean[np.newaxis,:])**2, weights=weights, axis=0 )
#    cvvar *= float(ndata_sets) / (ndata_sets - 1)
##    cvstd = np.std(heat_capacities, axis=0)
#    return cvmean, np.sqrt(cvvar)
        
    



def makeCvPlots(NDOF, eqTimeStep, basedirlist=['.'], ncurves=5, useWHAM=False,
                estimate_errors=True, neglect=0.001, verbose=False):
    """make a collection of Cv plots"""
    outprefix = "Cv"

    # load the data from the Visits.his files
    Vcollection = utils.VisitsHisCollection(basedirlist)
    Vcollection.check_data()

    Nrep = Vcollection.get_nreplicas() 
    Nvisits = Vcollection.get_nvisits()
    print "Number of replicas = ", Nrep 
    print "Number of histograms = ", Nvisits 

    # get the temperatures from file "temperatures" 
    Tlist = getTlist()

    # discard all VisitsHis objects with time before eqTimeStep
    Vcollection.tmax_discard(eqTimeStep)
    print "Using equilibration time step", Vcollection.vseries_list[0][0].time 
    
    Nvisits_init = Vcollection.get_nvisits()

    # reduce the data set to several equi-spaced time segments  
    indices = Vcollection.get_equispaced_indices(ncurves)
    Vcollection.reduce(indices)
    Vcollection.check_data()
    indices = list(reversed(range(Vcollection.get_nvisits())))

    Nvisits = len(indices)
    print "reducing", Nvisits_init, "datasets down to", Nvisits, "evenly distributed data sets"
    # indices[Nvisits-1] is the longest run

    if not useWHAM:
        neglect = None

    # estimate the error
    if estimate_errors:
        if Nvisits > 2:
            print "\nestimating errors in heat capacity"
            cvmean, cverr = estimate_error_jackknife(Vcollection, Tlist, NDOF)
        else:
            print "not enough datasets to estimate the error"
            estimate_errors = False

    #calculate all cv data
    print "\ncomputing heat capacity from", Nvisits, "Visits.his files"
    datalist = []
    Vall = Vcollection.vseries_list 
    ind = indices[0] # the longest run
    time1 = Vall[0][ind].time
    # main loop over indices data
    for isub in indices[1:]:
        # compute the cv for all data collected between time2 and time1
        time2 = Vall[0][isub].time
        # print "time 2 = ", time2
        print "computing Cv for the time period", time2, "to", time1 

        # visits corresponding to isub will be subtracted from ind  
        binenergy, visits = makeHistogram(Vall, ind, isub, neglect=neglect, verbose=verbose)
        
        assert(not np.any(np.isnan(visits)))
        if useWHAM:
            wham = WHAM.Wham1d(Tlist, binenergy, visits)
            wham.minimize()
            cvdata = wham.calc_Cv(NDOF, ntemp=500)
        else:
            cvdata = calcCvNoWHAM(binenergy, visits, Tlist, NDOF)
        assert(not np.any(np.isnan(cvdata)))
        #checkWHAM(binenergy, visits, cvdata, Tlist)
        datalist.append( (time1, time2, cvdata) )

    if len(indices) < ncurves:
        # plot Cv for the last Visits file without any subtraction.
        # Set time2 to None to indicate that no subtraction was done
        print "computing Cv from all existing data" 
        binenergy, visits = makeHistogram(Vall, ind, None)
        if useWHAM:
            wham = WHAM.Wham1d(Tlist, binenergy, visits)
            wham.minimize()
            cvdata = wham.calc_Cv(NDOF, ntemp=500)
        else:
            cvdata = calcCvNoWHAM(binenergy, visits, Tlist, NDOF)
        datalist.append( (time1, None, cvdata) )


    if True:
        # save the data to a file
        for time1, time2, cvdata in datalist:
            if time2 is None:
                fname = "%s%d" % (outprefix, time1)
            else:
                fname = "%s%d-%d" % (outprefix, time1, time2)
            print 'writing Cv data to file', fname
            np.savetxt(fname, cvdata)

        # save the data from the error estimate separately
        # because the temperatures are probably different            
        if estimate_errors:
            fnameerr = outprefix + ".err"
            print "writing error estimates to file", fnameerr
            with open(fnameerr, "w") as fout:
                fout.write("#T cv error\n")
                for data in itertools.izip(Tlist, cvmean, cverr):
                    fout.write("%g %g %g\n" % data)
                

    # now plot data
    fname=outprefix
    print 'saving Cv plot to file', fname
    for t1, t2, cvdata in datalist:
        if t2 is None:
            label = "%.2g" % (t1)
        else:
            label = "%.2g-%.2g=%.2g" % (t1, t2, t1-t2)
        plt.plot(cvdata[:,0], cvdata[:,5], "-", label = label)
    
    if estimate_errors:
        plt.errorbar(Tlist, cvmean, cverr, fmt=None, label="jackknife: direct (no wham)")
    plt.xlabel("T")
    plt.ylabel("Cv")
    leg = plt.legend(loc='best', fancybox='True')
    leg.get_frame().set_alpha(0.5)
    fname += ".pdf"
    with open(fname,"w") as fout:
        plt.savefig(fout, format="pdf")


#def usage():
#    print sys.argv[0], " -n natoms [-e equiTime] [dir1 dir2 ...]"
#    print " ------------------------------------------------------------------------------ " 
#    print " Make Cv plots from Visits.his.###### files and save to file Cv.pdf"
#    print "  The Visits.his.#### files are read from base directories dir1, dir2, etc."
#    print "  If no directories are supplied dir1 will be set to '.'"
#    print "  In each base directory the Visits.his.###### files are in subdirectories 1/ 2/ 3/ ..."
#    print "  The errors are estimated at each of the replica temperatures using the bootstrapping method"
#    print "  If many Visits.his.### files exists, only 5, evenly distributed, will be used"
#    print "    -n ndof:      number of degrees of freedom, required argument"
#    print "    -e timestep:  Visits.his.timestep would be subtracted from following ones to account for equilibration" 
#    print "                     optional, defaults to first Visits.his file " 
#    print "    -N:           Don't use WHAM. Calculate Cv only for the simulation temperatures."
#    print "    -b:           don't compute the estimated errrors"
#    print " ------------------------------------------------------------------------------ " 


def main():
    ############################################################
    # read command line parameters
    ############################################################
    import argparse
    parser = argparse.ArgumentParser(
        description="Make heat capacity plots from #/Visits.his.###### files and save to file Cv.pdf.  "
          "The histogram information is read from files basedir/[r]/Visits.his.[#] "
          "where basedir is '.' by default and [r] is in (1/, 2/, 3/, ... ) up to the number of replicas."
          "The errors are estimated at each of the replica temperatures using the bootstrapping resampling method."
          "If many Visits.his.### files exists, only 5, evenly distributed, will be used."
                                     )
    parser.add_argument("basedir", type=str, nargs="*", help="The simulation directories.  It is possible to specify"
                        " multiple base directories in order to accommodate resumed simulations."
                        "  The default is the working directory."
                        , default=["."])
    parser.add_argument("-n", "--ndof", type=int, help="Number of degrees of freedom.  "
                        "Supplying the number of degrees of freedom simply applies a constant shift to the heat capacity.  "
                        "This accounts for the momentum degrees of freedom"
                        " which are not explicitly simulated in a Monte Carlo simulation."
                        , default=0)
    parser.add_argument("-e", "--eq-timestep", type=int, help="Data before this time step will be ignored"
                        , default=-1)
    parser.add_argument("--no-wham", action="store_true", 
                        help="Don't use WHAM. Calculate Cv only for the simulation temperatures.")
    parser.add_argument("--no-errors", action="store_true", 
                        help="Don't compute the estimated errors")
    parser.add_argument("--neglect", type=float, default=0.01,
                        help="neglect bins which have fewer visits than neglect*max_visits, where max_visits is the maximum visits for that replica")
    parser.add_argument("--verbose", action="store_true", 
                        help="turn on verbose status printing")

    args = parser.parse_args()
    print args
    
    for d in args.basedir:
        assert os.path.isdir(d)
        
    dirlist = args.basedir
    for i in range(len(dirlist)):
        if dirlist[i].endswith('/'):
            dirlist[i] = dirlist[i][:-1]

    if args.neglect <= 0:
        args.neglect = None

    estimate_errors = not args.no_errors
    useWHAM = not args.no_wham    
    makeCvPlots(args.ndof, args.eq_timestep, basedirlist=dirlist, useWHAM=useWHAM, 
                estimate_errors=estimate_errors, neglect=args.neglect, verbose=args.verbose)

if __name__ == "__main__":
    main()
