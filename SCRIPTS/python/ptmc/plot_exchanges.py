import matplotlib as mpl
mpl.use("Agg") # so we can use it without an X server
import matplotlib.pyplot as plt

import sys
import os
import copy
import argparse
import numpy as np
import pickle
import multiprocessing as mp
from matplotlib.backends.backend_pdf import PdfPages

import ptmc_utils as utils 

from _exchanges_utils import ReadPermutations

def getReplicaPath(fname = "GMIN_out.1", dirlist=["."], makeplots=True):
    """
    get the paths of the replicas from file fname

    the replicas are undergoing permutations and are moving around between the temperature
    bins.  some definitions::

        position      : the index of the temperature bin
        replica label : the label of a replica

    Returns
    -------
    permutations : array, shape(nexchanges, nreplicas)
        permutations[t, i] is the previous position of the replica curently at
        position i at time t
    paths : array, shape(nexchanges, nreplicas)
        paths[t, r] is the position of replica r at time t
    replicas : array, shape(nexchanges, nreplicas)
        replicas[t, i] is the replica at position i at time t
    """
    temps = utils.getTemperatures()
    nreps = len(temps)
    readclass = ReadPermutations( nreps, dirlist=dirlist )
    readclass.readAll()
    #permutations, times = getPermutations(nreps, fname)
    permutations = np.array(readclass.permutations)
    times = np.array(readclass.times)
    nperm = len(permutations[:,0])
    print "nperm", nperm
    if True:
        #get `paths` and `replicas` from `permutations`
        try:
            from src.exch_tools import path_from_permutations
            paths, replicas, exchanges_counts = path_from_permutations(permutations)
        except ImportError:
            exchanges_counts = np.zeros(nreps, np.integer)
            opos = np.array(range(nreps))
            paths = np.zeros(np.shape(permutations), np.integer)
            replicas = np.zeros(np.shape(permutations), np.integer)
            paths[0,:] = range(nreps)
            replicas[0,:] = range(nreps)
            print "using slow python version of get paths"
            for i in range(1,nperm):
                perm = permutations[i-1,:]
                replicas[i,:] = replicas[i-1, perm]
                for pos,rep in enumerate(replicas[i,:]):
                    paths[i,rep] = pos
                exchanges_counts += (perm != opos)

    if True:
        """
        find how many round trips each replica made
        """
        print "calculating round trips"
        ndownup = np.zeros(nreps, np.integer)
        nupdown = np.zeros(nreps, np.integer)
        last_visit = np.zeros( nreps, np.integer)
        up = 1
        down = -1
        for i in range(nperm):
            iup = replicas[i,0]
            if last_visit[iup] == down:
                ndownup[iup] += 1
            last_visit[iup] = up

            idown = replicas[i,-1]
            if last_visit[idown] == up:
                nupdown[idown] += 1
            last_visit[idown] = down
        with open("exchanges.round_trips", "w") as fout:
            fout.write("#complete trips up and down\n")
            fout.write("#replica n_down_up n_up_down\n")
            for rep in range(nreps):
                print "complete trips up, down for replica %3d: %4d %4d" % ( rep, ndownup[rep], nupdown[rep])
                fout.write("%d %d %d\n" % (rep, ndownup[rep], nupdown[rep]))

    if not makeplots:
        print "not making plots"

    print "exchange frequency"
    ex_freq = exchanges_counts.astype(float) / nperm
    #permutations = np.array(permutations)
    #paths = np.array(paths)
    print exchanges_counts
    print ex_freq
    
    
    nperm = len(permutations[:,0])
    if False:
        with open("exchanges.perms", "w") as fout:
            print "writing file exchanges.perms"
            for i in range(nperm):
                ostr = ' '.join( [str(n) for n in permutations[i,:] ])
                fout.write(ostr + "\n")
    print "writing file exchanges.paths"
    np.savetxt("exchanges.paths", paths, fmt="%d")
        
    return permutations, paths, replicas, times, ex_freq

def makeplot_exchanges(paths, times, ex_freq, temps):
    print "plotting exchange paths"
    from matplotlib.collections import LineCollection
    nreps = len(temps)
    #paths = np.transpose(paths)
    print np.shape(paths)
    nexchanges = len(paths[:,0])
    timemin = times[0]
    timemax = times[-1]
    pp=PdfPages("exchanges.pdf")
    plt.subplot(1,2,1)
    plt.plot(ex_freq)
    plt.subplot(1,2,2)
    plt.plot(temps, ex_freq)
    pp.savefig()
    for i in range(nreps):
        plt.clf()
        temppath = np.array( [ temps[int(paths[t,i])] for t in range(nexchanges)])
        plt.plot( times, temppath, rasterized=True)
        # add a horizontal line at each temperature
        for T in temps:
            plt.plot( [timemin, timemax], [T,T], '-k')
        pp.savefig()
    pp.close()
    print "finished plotting exchange paths"

def makeplot_hist(permutations, temps):
    print "plotting exchange histograms"
    nreps = len(temps)
    print np.shape(permutations)
    pp=PdfPages("exchanges.hist.pdf")
    binedges = np.array([-0.5+i for i in range(nreps+1) ])
    binedgest = np.array([min(temps)] + [ np.mean(temps[i:i+1]) for i in range(nreps-1) ] + [max(temps)])
    xrange = [0,nreps]
    xranget = [min(binedgest), max(binedgest)]
    widths = abs( binedgest[1:] - binedgest[:-1] ) * 0.90
    for i in range(nreps):
        hist, b = np.histogram(permutations[:,i], binedges)
        norm = np.max(hist)
        plt.clf()
        plt.bar(binedgest[:-1], hist.astype(float)/norm, width = widths)
        #plt.plot(binedgest[:-1], hist.astype(float)/norm)
        plt.xlim(xranget)
        pp.savefig()
    pp.close()
    print "finished plotting exchange histograms"


def calculate_occupation_probability(paths, makeplots=True):
    """
    Parameters
    ----------
    paths : 2d array, shape(nexchanges, nreplicas)
        paths[t, r] is the position of replica r at time t
    
    Returns
    -------
    occupation_prob : 2d array, shape(nreplicas, nreplicas)
        occupation_prob[r, i] is the amount of time replica r spent at position i
    """
    print "calculating and plotting occupation probabilities"
    nexchanges, nreplicas = paths.shape

    try:
        from src.exch_tools import calculate_occupation_probability_cython
        print "using cython occupation probability method"
        occupation_prob = calculate_occupation_probability_cython(paths)
    except ImportError:
        print "using slow python occupation probability method"
        occupation_prob = np.zeros([nreplicas, nreplicas])
        for path in paths:
            for r, i in enumerate(path):
                occupation_prob[r, i] += 1.
    occupation_prob /= (float(nexchanges) / nreplicas) 
    
    # compute the "entropy" of each replica.  this is a measure of how flat the curves are
    flatness = np.nansum(np.log(occupation_prob), axis=1)
    flatness_perfect = 0.#nreplicas * np.log(1./nreplicas)

    if makeplots:
        pp=PdfPages("exchanges.occupation_prob.pdf")
        pmax = np.nanmax(np.nanmax(occupation_prob))
        for r, prob in enumerate(occupation_prob):
            plt.clf()
            plt.ylabel("occupation probability for replica %d" % r)
            plt.xlabel("temperature index (position)")
            plt.title("flatness = %g, should be %g" % (flatness[r], flatness_perfect))
            plt.plot(prob)
            plt.xlim([0,nreplicas-1])
            #plt.ylim([0,1])
            pp.savefig()
        pp.close()
    print "calculating and plotting occupation probabilities"


    

def calculateDiffusion(paths, makeplots=True):
    """
    calculate some sort of diffusion metric

    Parameters
    ----------
    paths : 2d array
        paths[t, r] is the position of replica r at time t

    Notes

    -----
    diffusion[p1, p2] : the average time elapsed since the replica at position
        p1 has been at p2.  time is measured in number of exchanges.

    rdiff[r, p] : the time elapsed since replica r was at position p
    """
    print "calculating and plotting diffusion of the replicas"
    #paths = paths.copy()
    #paths = paths[:1000,:]
    nexchanges, nreplicas = paths.shape
    diffusion = np.zeros([nreplicas, nreplicas])
    rdiff = np.ones([nreplicas, nreplicas], np.integer)
    maxelapsed_cum = 0
    maxelapsed = rdiff[0,0]
    with open("test.pkl", "wb") as fout:
        pickle.dump(diffusion, fout)

    try:
#        raise ImportError()
        from src.exch_tools import calculate_diffusion_cython
        print "calculating diffusion using c module"
        diffusion1, rdiff, maxelapsed, maxelapsed_cum = calculate_diffusion_cython(paths)
        
    except ImportError:
        print "calculating diffusion using slow python code"
        for t, positions in enumerate(paths):
            maxelapsed += 1
            maxelapsed_cum += maxelapsed
            rdiff[:,:] += 1
            for r, p in enumerate(positions):
                rdiff[r, p] = 0
    
            for r in xrange(nreplicas):
                p1 = paths[t, r]
                for p2, telapsed in enumerate(rdiff[r,:]):
                    diffusion[p1, p2] += telapsed

    # deal with situations where replica r was never at position j
    # set these to np.nan.  The array must be floats for this to work
    diffusion = np.array(diffusion1, dtype=np.float64)
    indices = np.where(abs(diffusion - maxelapsed_cum) < 1e-5)
    diffusion[indices] = np.nan
    diffusion /= nexchanges
        

    if False: #debugging stuff
        print ""
        print "rdiff"
        for rd in rdiff:
            print rd

        print ""
        print "diffusion" 
        print ""
        for rd in diffusion:
            print rd
        print maxelapsed, maxelapsed_cum
        print np.max(diffusion) * nexchanges
        print nreplicas * nreplicas


    with open("test.pkl", "wb") as fout:
        pickle.dump(diffusion, fout)

    if makeplots:
        pp=PdfPages("exchanges.diffusion.pdf")
        ymax = np.nanmax(diffusion)
        for p, d in enumerate(diffusion):
            plt.clf()
            plt.ylabel("average time to get to %d (time in number of exchanges)" % p)
            plt.xlabel("temperature index")
            plt.plot(d)
            plt.xlim([0,nreplicas-1])
            plt.ylim([0,ymax])
            pp.savefig()
        pp.close()
    print "finished calculating and plotting diffusion of the replicas"

def main():
    ############################################################
    # read command line parameters
    ############################################################
    parser = argparse.ArgumentParser(
        description="Get information about the PT exchanges and make plots showing replica paths through "
                    "temperature space, replica diffusion, and exchange probabilities, and occupation probabilities."
                    "The exchange information is read from file exchanges.permutations.  If this file does not exist"
                    "I try to get the exchange information from GMIN_out.1.  "
                    "Some of these routines can be quite slow when done in pure python, so there is an extension module"
                    "written in c.  This must to be compiled.  "
                    "Furthermore, by default the analysis is done in parallel with 4 cores.  This can be turned off with the options"
                    )
    parser.add_argument("basedir", type=str, nargs="*", help="The simulation directories.  It is possible to specify"
                        " multiple base directories in order to accommodate resumed simulations."
                        "  The default is the working directory."
                        , default=["."])
    parser.add_argument("--serial", action="store_true", help="Do the computations in serial, not in parallel.")
    args = parser.parse_args()
    print args
    
    for d in args.basedir:
        assert os.path.isdir(d)
        
    dirlist = args.basedir
    for i in range(len(dirlist)):
        if dirlist[i].endswith('/'):
            dirlist[i] = dirlist[i][:-1]


    temps = utils.getTemperatures(dirlist[0])

    # read the permutations from disc and change the formatting into replica paths.
    permutations, paths, replicas, times, ex_freq = getReplicaPath(dirlist=dirlist, makeplots=True)
    
    # set up the parallel computation
    if not args.serial:
        ncores = 4
        print "making a pool of workers"
        pool = mp.Pool(processes=ncores)
        myapply = pool.apply_async
    else:
        myapply = apply
        ncores = 1
    
    # make some plots and do some calculations on the replica paths
    myapply(makeplot_exchanges, (paths, times, ex_freq, temps))
    
    myapply(makeplot_hist, (permutations, temps))
    myapply(calculate_occupation_probability, (paths,))
    myapply(calculateDiffusion, (paths,))

    if ncores > 1:
        print "finished submitting all the jobs in parallel"
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()
