import itertools
import os
import gzip

#from numpy import *
import numpy as np



class ReadDataFile:
    def __init__(self, datafile, natoms=None):
#        self.frozen=zeros(natoms, integer)
        self.freeze = False
        self.frozenlist = []
        sp0 = []
        with open(datafile,"r") as f:
            for line in f:
                sp = sp0 + line.split()

                if sp[-1] == "+++": # +++ means the line is continued on the next line 
                    sp0 = sp[:-1]
                    continue
                else:
                    sp0 = []                
                
                if len(sp) == 0: continue
                if sp[0] == 'FREEZE':
                    self.freeze = True
                    for i in sp[1:]:
                        self.frozenlist.append(int(i) - 1)
#                        self.frozen[int(i)-1] = 1
                elif sp[0] == 'PERIODIC':
                    self.boxl = float(sp[1])
                elif sp[0] == 'BINARY':
                    self.ntypeA = int(sp[1])
                    self.epsAB = float(sp[2])
                    self.epsBB = float(sp[3])
                    self.sigAB = float(sp[4])
                    self.sigBB = float(sp[5])
                elif sp[0] == 'SHIFTCUT':
                    self.rcut = float(sp[1])
                elif sp[0] == 'RESTRICTREGION':
                    self.restrict_region_radius = float(sp[1])
                elif sp[0] == 'BSPT':
                    self.histmin = float(sp[1]) #  quench histogram
                    self.histmax = float(sp[2]) #  quench histogram
                    self.ptemin = float(sp[3]) #  instantaneous energy histogram
                    self.ptemax = float(sp[4]) #  instantaneous energy histogram
                    self.nenrper = int(sp[11]) #  instantaneous energy histogram
                    self.hbins = int(sp[12]) #  quench histogram
        
                elif sp[0] == 'PTMC':
                    self.histmin = float(sp[1]) #  quench histogram
                    self.histmax = float(sp[2]) #  quench histogram
                    self.ptemin = float(sp[3]) #  instantaneous energy histogram
                    self.ptemax = float(sp[4]) #  instantaneous energy histogram
                    self.exrate = float(sp[7].replace("D","E"))
#                    self.nenrper = int(sp[11]) #  instantaneous energy histogram
#                    self.hbins = int(sp[12]) #  quench histogram

        
        if self.freeze:
            self.frozenlist = np.array(self.frozenlist)
            self.nfrozen = len(self.frozenlist)
            
            if natoms is not None:
                self.natoms = natoms
                self.frozen = np.zeros(natoms, np.integer)
                self.frozen[self.frozenlist] = 1
                
                frozen = self.frozen

                self.nmobile = natoms-self.nfrozen
                self.mobilelist = np.where(frozen[:]==0)[0]
                self.frozenlist = np.where(frozen[:]!=0)[0]
        
                if hasattr(self, "ntypeA"):        
                    ntypeA = self.ntypeA
                    self.mobileA =    np.where(frozen[0:ntypeA]==0)[0]
                    self.mobileB =    np.where(frozen[ntypeA:]==0)[0]
                    self.mobileB += ntypeA
                    self.frozenA =    np.where(frozen[0:ntypeA]!=0)[0]
                    self.frozenB =    np.array([i for i in self.frozenlist if i >= ntypeA], np.integer)

class read_data_file(ReadDataFile):
    """for backward compatibility"""
    pass


def read_Visits(fname):
    """read a Visits.his file"""
    elist = []
    countlist = []
    if fname.endswith('.gz'):
        fin = gzip.open(fname, "r")
    else:
        fin = open(fname, "r")
    #with myopen(fname, "r") as fin:
    if True:
        #skip until the first word is Visits
        reading = False
        for line in fin:
            sline = line.split()
            if reading:
                if sline[0] == "Visits":
                    break
                elist.append( float(sline[0]) )
                countlist.append( float(sline[1]) )

            if sline[0] == "Visits":
                reading = True
    fin.close()
    return np.array(elist), np.array(countlist)

def getTemperatures(directory="."):
    temps = np.loadtxt(directory + "/temperatures")
    return temps

def get_dir_list(basedir='.'):
    """return list of directories 1/ 2/ 3/ ...
    """
    dirlist=[]
    started = False
    for i in range(100):
        basename = str(i)
        dir0 = basedir + '/' + basename
        if os.path.isdir(dir0):
            dirlist.append(basename)
            started = True
        elif started:
            return dirlist
    raise Exception("error: didn't find any directories")
    return dirlist


class VisitsHis(object):
    """an object to represent a single Visits.his file
    """
    def __init__(self, fname, dir="."):
        self.fname = os.path.basename(fname)
        if self.fname != fname:
            self.dir = os.path.dirname(fname)
            if dir != ".":
                print "VisitsHis: warning: ", fname, ": setting dir to", self.dir, ". Ignoring", dir
        else:
            self.dir = dir
        self.fullfname = self.dir+"/"+self.fname
        #get the time from the file name 
        #  Visits.his   the final output file
        #  Visits.his.######.0   the #### gives the time
        #  first strip of the .gz if it exists
        newfname = self.fname
        if newfname.endswith( '.gz'):
            newfname = newfname[:-3]
        if newfname == "Visits.his":
            #the final file
            self.time = None
        elif newfname == "Visits.his.restart":
            #the final file
            self.time = None
        else:
            self.time = int(float( newfname[ len("Visits.his.") : ] ))

    def read(self):
        if not hasattr(self, "counts_all"): 
            self.energies_all, self.counts_all = read_Visits(self.fullfname)
            self.nonzero = np.where(self.counts_all[:] != 0)[0]
            self.energies = self.energies_all[self.nonzero]
            self.counts = self.counts_all[self.nonzero]

    def getHis(self):
        return self.energies, self.counts

    def getHisAll(self):
        return self.energies_all, self.counts_all

    def getHisNorm(self):
        return self.energies, self.counts / sum(self.counts)

    def getHisDiffAll(self, rep2):
        c = self.counts_all - rep2.counts_all
        return self.energies_all, c

    def getHisDiff(self, rep2):
        #c = self.counts_all - rep2.counts_all
        e, c = self.getHisDiffAll( rep2 )
        nz = np.where(c[:] != 0)[0]
        return e[nz], c[nz]

    def getHisDiffNorm(self, rep2):
        e, c = self.getHisDiff(rep2)
        return e, c / sum(c)

class VisitsHisSeries(list):
    """a list of VisitsHis files generated at successive times from the same temperature
    """
#    def __init__(self):
#        self.vlist = []
    
    def mysort(self):
        self.sort(key = lambda v:v.time)
        
    
    def read_directory(self, d):
        """
        read Visits.his.* files in a directory and add them to the total list
        """
        files = os.listdir(d)
        f = "Visits.his."
        pre = 'Visits.his.'
        newvlist = []
        for s in files:
            if s.startswith(pre):
                if s == 'Visits.his.gz': continue
                v = VisitsHis(s, d)
                if v.time is not None:
                    newvlist.append(v)
        self += newvlist
        self.mysort()
    
    def get_nvisits(self):
        return len(self)
    
#    def __getitem__(self, i):
#        return self.vlist[i]

class VisitsHisCollection(object):
    """
    organize all the Visits.his files for a parallel tempering run
    """
    def __init__(self, basedirlist=["."]):
        self.basedirlist = basedirlist
        self.dirlist = get_dir_list()
        self.vseries_list = []
        
        self.get_all_visits(basedirlist)
    
    def get_all_visits(self, basedirlist):
        for d in self.dirlist:
            vseries = VisitsHisSeries()
            self.vseries_list.append(vseries)
            
            for basedir in basedirlist:
                vseries.read_directory(basedir + "/" + d)
    
    def remove_duplicates(self):
        for r in range(self.get_nreplicas()):
            # eliminate the duplicates
            vseries_dict = dict([(v.time, v) for v in self.vseries_list[r]])
            # rebuild the VisitsHisSeries
            self.vseries_list[r] = VisitsHisSeries(vseries_dict.values())
            self.vseries_list[r].mysort()

    
    def remove_uncommon_times(self):
        timeset = set()
        for vs in self.vseries_list:
            times = [v.time for v in vs]
            timeset.update()
    
    def check_data(self):
        self.remove_duplicates()
        
        nvisits = self.get_nvisits()
        for vseries in self.vseries_list[:]:
            if nvisits != vseries.get_nvisits():
                raise Exception("number of VisitsHis files are not equal %d %d" % (nvisits, vseries.get_nvisits()))
        
        
        vs0 = self.vseries_list[0]
        for vs in self.vseries_list[1:]:
            for v, v0 in itertools.izip(vs, vs0):
                if v.time != v0.time:
#                    print "VisitsHis objects have different times %s %s", (str(v.time), str(v0.time))
#                    self.remove_uncommon_times()
#                    self.check_data()
#                    return
                    raise Exception("VisitsHis objects have different times %s %s", (str(v.time), str(v0.time)))
        
        print "VisitsHisCollection: check ok"
    
    def tmax_discard(self, tmax):
        """
        discard all data that has time < tmax
        """
        for i in range(len(self.vseries_list)):
            vs = self.vseries_list[i]
            self.vseries_list[i] = VisitsHisSeries([v for v in vs if v.time > tmax])
#            imax = 
            
#    def get_Vall(self):
#        return self.vseries_list
    def get_nreplicas(self):
        return len(self.vseries_list)
    def get_nvisits(self):
        return self.vseries_list[0].get_nvisits()
    def get_time(self, index):
        return self.vseries_list[0][index].time
        
    
    def get_equispaced_indices(self, nnew=5):
        """return a set roughly nnew of indices roughly equispaced in time
        
        the spacing is not guaranteed to be uniform
        
        the number of indices returned will be as close as possible to nnew.  It will
        be between 0 and 2*nnew 
        """
        ntot = self.get_nvisits()
        if ntot < 2*nnew:
            #indices = reversed(range(ntot))
            indices = range(ntot-1,-1,-1)
        else:
            dn = ntot / nnew
            indices = range(ntot-1, -1, -dn)
        return list(reversed(indices))
    
    def reduce(self, indices):
        """
        remove all data from the visits series except those given by indices
        """
        for j in range(len(self.vseries_list)):
            vs = self.vseries_list[j]
            self.vseries_list[j] = VisitsHisSeries([vs[i] for i in indices])


        
            
        

        

class SystemHis(object):
    """
    a collection of VisitsHis objects which define the state of the system at a given time
    """
    def __init__(self, dirs, time):
        self.dirs = dirs
        self.time = time
        self.getVisits(dirs, time)

    def getVisits(self, dirs, time):
        self.visitslist = []
        for d in dirs:
            fname = "Visits.his." + str(time) + ".0"
            self.visitslist.append( VisitsHis( fname, d ) )


def getNReplicas():
    """
    try to figure out the number of replicas
    """
    nreps = None
    #count the number of lines in file temperatures
    try:
        with open("temperatures", "r") as fin:
            nreps = 0
            for line in fin:
                nreps += 1
        return nreps
    except:
        pass
    #count directories named 1/, 2/, 3/, ... 
    if True:
        nreps = 1
        while os.path.isdir(str(nreps)):
            nreps += 1
        nreps -= 1
        if nreps > 1:
            return nreps
        else:
            nreps = None
    print "warning: couldn't determine the number of replicas"
    return None


class Visits2His(object):
    """
    read the binary data in Visits2.his files and convert it into a 2d array
    
    Parameters
    ----------
    fname : str
        the name of the file to read.  Can be a name or a full path
    dir : str, optional
        the directory that the file is in.  this option is only applicable if
        fname is simply a file name and not a path
    """
    def __init__(self, fname, dir="."):
        self.fname = os.path.basename(fname)
        if self.fname != fname:
            self.dir = os.path.dirname(fname)
            if dir != ".":
                print "VisitsHis: warning: ", fname, ": setting dir to", self.dir, ". Ignoring", dir
        else:
            self.dir = dir
        self.fullfname = self.dir+"/"+self.fname
        #get the time from the file name 
        #  Visits.his   the final output file
        #  Visits.his.######.0   the #### gives the time
        #  first strip of the .gz if it exists
        newfname = self.fname
        if newfname.endswith( '.gz'):
            raise Exception("Visits2His can't deal with zipped files yet")
            newfname = newfname[:-3]
        if newfname == "Visits2.his":
            #the final file
            self.time = None
        elif newfname == "Visits2.his.restart":
            #the final file
            self.time = None
        else:
            self.time = int(float( newfname[ len("Visits2.his.") : ] ))

    def read(self):
        """
        read the data from the binary file dumped with fortran WRITE 
        with format UNFORMATTED
        
        Returns
        -------
        visits : 1d array
            this is the flattened 2d visits array.  it can be reshaped to have
            shape (nqbins, nibins) where nqbins is the number of quenched energy bins and
            nibins is the number of instantaneous energy bins
        """
        from fortranfile import FortranFile
        f = FortranFile(self.fullfname)
        data = f.readInts()
        self.visits = data
        return data
        

def wstd_np(x, w):
    mean, wsum = np.average(x, weights=w, returned=True)
    mean2 = np.average(x**2, weights=w)
    var = (mean2 - mean**2) #/ wsum
#    cv2mean = np.average(heat_capacities**2, weights=weights, axis=0)
    return mean, np.sqrt(var)

def wstd_average(x, w):
    """return the weighted mean and the standard deviation of the weighted mean"""
    t = np.sum(w)
    t2 = np.sum(w**2)
    xbar = np.sum(w * x) / t
    dx2 = np.sum(w * (x-xbar)**2)
    sig2 = t/(t**2 - t2) * dx2
    return xbar, np.sqrt(sig2)

def wstd(x, w):
    """
    return weighted mean and standard deviation of the samples
    """
    mean, wsum = np.average(x, weights=w, returned=True)
    mean2 = np.average(x**2, weights=w)
    var = (mean2 - mean**2) #/ wsum
#    cv2mean = np.average(heat_capacities**2, weights=weights, axis=0)
    return mean, np.sqrt(var)



