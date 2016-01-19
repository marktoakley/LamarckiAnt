import matplotlib as mpl
mpl.use("Agg") # so we can use it without an X server
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import os
import getopt, sys
import copy
from itertools import izip
import numpy as np

import ptmc_utils as utils

def listVisitsFiles(dir1):
    """
    return a sorted list of VisitsHis objects from files of the form
    Visits.his.#######
    """
    files = os.listdir(dir1)
    f = "Visits.his."
    #files = filter( lambda s: ( f in s ), files)
    pre = 'Visits.his.'
    visits = []
    for s in files:
        if pre in s:
            #print s
            visits.append( utils.VisitsHis( s , dir1) )
    visits = sorted(visits, key=lambda a: a.time)
    return visits

def reduce_list(nums, nnew=5):
    nfiles = len(nums)
    if nnew != None:
        if nfiles < 2*nnew: return nums
        frq = nfiles/nnew
        nums = [ nums[nfiles-i-1] for i in reversed(range(0,len(nums),frq)) ]
    #print nums
    return nums

def reduce_list_indices(nums, nnew=5):
    nfiles = len(nums)
    indices = reversed(range(nfiles))
    if nnew != None:
        if nfiles > 2*nnew: 
            frq = nfiles/nnew
            #nums = [ nums[nfiles-i-1] for i in reversed(range(0,len(nums),frq)) ]
            indices = [ nfiles-i-1 for i in reversed(range(0,nfiles,frq)) ]
    return indices

def getDirList( basedir = "." ):
    dirlist=[]
    started = False
    for i in range(100):
        dir = "%s/%d" % (basedir, i)
        if os.path.isdir(dir):
            dirlist.append(dir)
            started = True
        elif started:
            return dirlist
    print "error: didn't find any directories"
    return dirlist

def plotVisits(nfiles=5, basedirlist = ['.']):

    if False:
        dirlist = getDirList()
        print dirlist

        files = listVisitsFiles(dirlist[0])
        files = reduce_list(files, nfiles)

        file = files[-1]

        allvisits = []
        for dir in dirlist:
            allvisits.append( listVisitsFiles( dir ) )

    #get all visits files
    Vcollection = utils.VisitsHisCollection(basedirlist)
    Vcollection.check_data()

    #we don't want to plot them all, so select nfiles objects
    indices = Vcollection.get_equispaced_indices(nfiles)
    Vcollection.reduce(indices)
    Vcollection.check_data()
    indices = list(reversed(range(Vcollection.get_nvisits())))
    allvisits = Vcollection.vseries_list

    print "indices", indices

    if True:
        pp=PdfPages("Visits.his.pdf")

        plt.xlabel("energy")
        plt.ylabel("visits")
        for visits in allvisits:
            myvis = visits[-1]
            #n,fname = file
            #myvis = utils.VisitsHis(fname, dir)
            myvis.read()
            e, counts = myvis.getHis()
            plt.plot((e), (counts), '-')
        pp.savefig()
        plt.clf()
        
        plt.xlabel("energy")
        plt.ylabel("log_10 visits")
        for visits in allvisits:
            myvis = visits[-1]
            #n,fname = file
            #myvis = utils.VisitsHis(fname, dir)
            myvis.read()
            e, counts = myvis.getHis()
            plt.plot((e), np.log10(counts), '-')
        pp.savefig()
        plt.clf()

        plt.xlabel("energy")
        plt.ylabel("visits")
        for visits in allvisits:
            myvis = visits[-1]
            oldvis = visits[ len(visits)/2 ]
            title = "%g - %g = %g"%(myvis.time, oldvis.time, myvis.time - oldvis.time)
            #n,fname = file
            #myvis = utils.VisitsHis(fname, dir)
            myvis.read()
            oldvis.read()
            e, counts = myvis.getHisDiff(oldvis)
            plt.plot((e), (counts), '-')
        plt.title(title)
        pp.savefig()
        plt.clf()
        
        plt.xlabel("energy")
        plt.ylabel("log_10 visits")
        for visits in allvisits:
            myvis = visits[-1]
            oldvis = visits[ len(visits)/2 ]
            title = "%g - %g = %g"%(myvis.time, oldvis.time, myvis.time - oldvis.time)
            #n,fname = file
            #myvis = utils.VisitsHis(fname, dir)
            myvis.read()
            oldvis.read()
            e, counts = myvis.getHisDiff(oldvis)
            plt.plot((e), np.log10(counts), '-')
        plt.title(title)
        pp.savefig()
        plt.clf()

        plt.ylabel("log density of states")
        plt.xlabel("log energy")
        Tlist = np.genfromtxt("temperatures")
        for T, visits in izip(Tlist, allvisits):
            myvis = visits[-1]
            oldvis = visits[ len(visits)/2 ]
            title = "%g - %g = %g"%(myvis.time, oldvis.time, myvis.time - oldvis.time)
            #n,fname = file
            #myvis = utils.VisitsHis(fname, dir)
            myvis.read()
            oldvis.read()
            e, counts = myvis.getHisDiff(oldvis)
            log_dos = np.log(counts) + e / T 
            plt.plot(e, log_dos, '-')
        pp.savefig()
        plt.clf()
        
        
        pp.close()




    Tlist = utils.getTemperatures()

    pp=PdfPages("Visits.his.all.pdf")
    plt.clf()
    rep = 0
    for T, visits in izip(Tlist, allvisits):
        rep += 1
        visits = reduce_list( visits, nfiles)
        finalvis = visits[-1]
        finalvis.read()
        for myvis in reversed(visits[:-1]):
            title = "replica %d, T = %f" % (rep, T)
            plt.title( title)

            myvis.read()

            #e, counts = utils.read_Visits(fname)
            e, counts = finalvis.getHisDiffNorm(myvis)
            label = "%g - %g = %g"%(finalvis.time, myvis.time, finalvis.time - myvis.time)
            label = "%g"%(finalvis.time - myvis.time)
            plt.ylabel("log_{10} Num. Visits")
            plt.plot(e, np.log10(counts), '-', label=label)
        plt.legend(loc="lower left")
        pp.savefig()
        plt.clf()
    pp.close()


def usage():
    print sys.argv[0], " [-n nfiles] [dir1 dir2 ...]"
    print "Make visits plots from Visits.his.###### files and save to file Visits.his.pdf"
    print "  The Visits.his.#### files are read from base directories dir1, dir2, etc."
    print "  If no directories are supplied dir1 will be set to '.'"
    print "  In each base directory the Visits.his.###### files are in subdirectories 1/ 2/ 3/ ..."
    print "    -n nfiles:  number of curves on each plot"

if __name__ == "__main__":
    ############################################################
    # read command line parameters
    ############################################################
    import argparse
    parser = argparse.ArgumentParser(
        description="Make plots of the visits histograms from #/Visits.his.###### files and save to file Visits.his.pdf.  "
          "The histogram information is read from files basedir/[r]/Visits.his.[#] "
          "where basedir is '.' by default and [r] is in (1/, 2/, 3/, ... ) up to the number of replicas."
          "If many Visits.his.### files exists, only 5, evenly distributed, will be used."
                                     )
    parser.add_argument("basedir", type=str, nargs="*", help="The simulation directories.  It is possible to specify"
                        " multiple base directories in order to accommodate resumed simulations."
                        "  The default is the working directory."
                        , default=["."])
    args = parser.parse_args()
    print args


    dirlist = args.basedir
    for d in dirlist:
        assert os.path.isdir(d)
    for i in range(len(dirlist)):
        if dirlist[i].endswith('/'):
            dirlist[i] = dirlist[i][:-1]


    nfiles = 5

    plotVisits(nfiles=nfiles, basedirlist=dirlist)
