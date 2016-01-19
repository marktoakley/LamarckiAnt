
import gzip
import os
import numpy as np


class ReadExchanges(object):
    """
    get exchanges by reading GMIN_out.## files
    
    this class is deprecated but still exists for backwards compatibility.  use ReadPermutations instead

    ensure ability to combine exchanges from multiple GMIN_out.1 files from
    runs that have been restarted

    This will be tricky because if a run is ended prematurely then many
    exchanges will be registered but then abandoned when the run is restarted
    from the last checkpoint.  I have to make sure to account for this correctly
    """
    def __init__(self, nreps, dirlist = ['.']):
        self.nreps = nreps
        self.dirlist = dirlist
        #remove trailing / from directories
        for i in range(len(dirlist)):
            if dirlist[i].endswith('/'):
                dirlist[i] = dirlist[i][:-1]
        #get times
        self.timelist = [ self.getEndTime(d) for d in self.dirlist ]
        print "times", self.timelist
        #sort the directories based on time
        timedir = sorted( zip(self.timelist, self.dirlist) )
        self.dirlist = [ z[1] for z in timedir ]
        self.timelist = [ z[0] for z in timedir ]
        print "timedir", timedir
        print "times", self.timelist

        self.exchanges = []


    def getEndTime(self, d):
        """
        try to read the file 1/bsptrestart to find the time when the simulation
        ended.  The time should be the first field in the first line of bsptrestart
        """
        with open( d+"/1/bsptrestart", "r") as fin:
            line = fin.readline()
            time = float(line.split()[0])
        return time

    def readAll(self):
        for time, d in zip(self.timelist, self.dirlist):
            print 'reading from directory', d
            exch = read_exchanges(self.nreps, d, time)
            self.exchanges += exch

        if True: #test
            print "checking exchanges are in sorted order"
            tprev = -1
            for ex in self.exchanges:
                if ex[0] <= tprev:
                    print "warning: exchanges are not in sorted order", tprev, ex[0]
                tprev = ex[0]
        return self.exchanges

def read_ex(fname, maxtime = None ):
    """parse the GMIN_out.## file to get the exchange information
    """
    try:
        fin = open(fname, "r")
    except:
        fin = gzip.open(fname+".gz", "r")
    exlist = []
    for line in fin:
        if "Exchange from this replica" in line:
            sline = line.split()
            rep1 = int(sline[5])
            rep2 = int(sline[10])
            e1 = float(sline[7])
            e2 = float(sline[12])
            time = float(sline[15])
            if maxtime != None:
                if time > maxtime: break
            exlist.append( (time, rep1, rep2, e1, e2) )
            #print time, rep1, rep2, e1, e2
    fin.close()
    return exlist

def read_ex_2( fname, maxtime = None ):
    """parse the GMIN_out.## file to get the exchange information

    some GMIN_out files have a different format
    time is at column 14 rather than 15
    """
    try:
        fin = open(fname, "r")
    except:
        fin = gzip.open(fname+".gz", "r")
    exlist = []
    for line in fin:
        if "Exchange from this replica" in line:
            sline = line.split()
            rep1 = int(sline[5])
            rep2 = int(sline[10])
            e1 = float(sline[7])
            e2 = float(sline[12])
            time = float(sline[14])
            if maxtime != None:
                if time > maxtime: break
            exlist.append( (time, rep1, rep2, e1, e2) )
            #print time, rep1, rep2, e1, e2
    fin.close()
    return exlist

def read_exchanges(nreps, d = '.', maxtime=None):
    """parse the GMIN_out.## file to get the exchange information
    
    this will call the appropriate function
    """
    exall = []
    for rep in range(nreps):
        fnamebase = d+"/GMIN_out." 
        fname = fnamebase + str(rep+1)
        try:
            exchanges = read_ex(fname, maxtime)
        except:
            exchanges = read_ex_2(fname, maxtime)
        exall += exchanges
    exall = sorted(exall) #sort by time
    return exall


class ReadPermutations(ReadExchanges):
    """
    get exchange permutations by reading file exchanges.permutations
    
    if that file is not present it will try to parse the required information from GMIN_out.1

    ensure ability to combine exchanges from multiple GMIN_out.1 files from
    runs that have been restarted

    This will be tricky because if a run is ended prematurely then many
    exchanges will be registered but then abandoned when the run is restarted
    from the last checkpoint.  I have to make sure to account for this correctly
    
    Parameters
    ----------
    nreps : integer
        number of replicas
    dirlist : list of strings
        list of directories that hold info about this run
    """
    def __init__(self, nreps, dirlist = ['.']):
        super(ReadPermutations, self).__init__(nreps, dirlist=dirlist)
        self.times = []
        self.permutations = []



    def readAll(self):
        for time, d in zip(self.timelist, self.dirlist):
            print 'reading from directory', d
            newperms, newtimes = getPermutations(self.nreps, directory = d)
            """ discard all permutations after time time """
            nperms = len(newperms)
            #print "readAll: nperms", nperms, time, newtimes[-5:]
            for i in reversed(range(nperms)):
                if newtimes[i] <= time:
                    #print "times:", i, newtimes[i],  time
                    imax = i+1
                    if imax < nperms:
                        newtimes = newtimes[:imax]
                        newperms = newperms[:imax]
                    break
            self.times += newtimes
            self.permutations += newperms
            #print len(self.times)

        if True: #test
            print "checking exchanges are in sorted order"
            tprev = -1
            for time in self.times:
                if time <= tprev:
                    print "warning: exchanges are not in sorted order", tprev#, ex[0]
                tprev = time

def getPermutations_new(directory):
    """try to read the exchanges from file exchanges.permutations
    """
    fname = directory + "/exchanges.permutations"
    fnamegz = fname + ".gz"
    myopen = open
    if not os.path.isfile(fname):
        if os.path.isfile(fnamegz):
            fname = fnamegz
            myopen = gzip.open
        else:
            raise IOError("file %s does not exist" % fname)

    
    print "reading exchanges from file", fname
    data = np.genfromtxt(myopen(fname, "r"))
    #data = np.genfromtxt(fname)
    times = data[:,0]
    permutations = data[:,1:] - 1 #  subtract 1 for fortran indexing

    return permutations, times
    

def getPermutations(nreps, fname="GMIN_out.1", directory="."):
    """
    return permutations and monte carlo steps at which they occur
    
    Returns
    --------
    permutations : list of integer numpy arrays
    times : list of numpy arrays
    """
    try:
        permutations, times = getPermutations_new(directory)
        lperms = [np.array(p, np.integer) for p in permutations]
        return lperms, list(times)
    except IOError:
        pass
    fname = directory + "/" + fname
    permutations = []
    times = []
    if True:
        try:
            fin = open(fname, "r")
        except:
            fin = gzip.open(fname+".gz", "r")
        for line in fin:
            sline = line.split()
            #try to get the time of the permutation
            if len(sline) >= 3:
                if sline[1] == "Vn=": #this will be used for old runs before the actual time was printed
                    time = float(sline[0])
                    continue
                if sline[0] == "replica" and sline[1] == "permutation":
                    time = float(sline[2])
                    continue
            #now get the permutation.  There's no unique way to identify the line, so try several things
            if len(sline) != nreps: continue
            try:
                i1 = int(sline[0])
            except ValueError:
                continue
            if i1 > nreps: continue
            #print line
            #time = int(sline[0])
            permutation = [ int(s) for s in sline[0:nreps] ]
            #permutation -= 1  #correct fortran indices
            permutations.append( np.array(permutation) - 1 )
            times.append(time)
        fin.close()
    return permutations, times
