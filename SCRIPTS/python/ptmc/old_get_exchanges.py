import gzip
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import ptmc_utils as utils
import plotVisits
import copy


#read PT exchanges from GMIN_out.## files

    
class RoundTrips(object):
    def __init__(self, rep, nreps):
        self.top = nreps-1
        self.bottom = 0
        self.ntrips_up = 0
        self.ntrips_down = 0
        self.status = rep
    def register(self, newrep):
        if newrep == self.top:
            if self.status == self.bottom:
                self.ntrips_up += 1
            self.status = self.top
        elif newrep == self.bottom:
            if self.status == self.top:
                self.ntrips_down += 1
            self.status = self.bottom
    def results(self):
        return self.ntrips_up, self.ntrips_down


def find_path(rep, exchanges, nreps):
    time0 = exchanges[0][0]
    path = [(time0, rep)]
    trip_record = RoundTrips( rep, nreps )
    for ex in exchanges:
        if ex[1] == rep:
            time = ex[0]
            rep = ex[2]
            trip_record.register(rep)
            path.append( (time,rep) )
        elif ex[2] == rep:
            time = ex[0]
            rep = ex[1]
            trip_record.register(rep)
            path.append( (time,rep) )
    return np.array(path), trip_record.results()

class MyEta():
    def __init__(self, temps, f):
        import scipy.interpolate
        self.finterp = scipy.interpolate.UnivariateSpline(temps, f, s=0)
        #self.fnorm = finterp.integral(temps[0], temps[-1])
        self.temps = temps
    def dfdT(self, T):
        #print "derivatives", T, self.finterp.derivatives(T)
        try:
            return np.array( [self.finterp.derivatives(t)[1] for t in T] )
        except:
            df = self.finterp.derivatives(T)[1]
            return df
    def _dT(self, T):
        """
        return the temperature interval around T,  return T_{i+1} - T_i
        T_i < T < T_{i+1}
        """
        if T == self.temps[0]: return self._dT(T+1e-15)
        if T == self.temps[-1]: return self._dT(T-1e-15)
        import bisect
        #returns the index of the first element large than T.  Will fail for Tmax 
        i = bisect.bisect(self.temps, T) 
        #print "T bisect", self.temps[i], T, self.temps[i-1]
        return self.temps[i] - self.temps[i-1]
    def dT(self, mytemps):
        try:
            return np.array( [self._dT(T) for T in mytemps] )
        except: 
            return self._dT(T)
    def eta(self, T):
        return np.sqrt( self.dfdT(T) / self.dT(T) )


def getDiffusion(exchanges, nreps, temps):
    """
    get information about the diffusion of replicas between the 
    highest and lowest temperatures
    
    D local diffusion rate
    j current
    eta = distribution of temperatures
    etaopt = optimal distribution of temperatures

    new temperatures are chosen from etaopt by
    integral_Tmin^Tk etaopt dT = k / nreps
    """
    nup, ndown = getNupNdown(exchanges, nreps)
    f = nup / (nup + ndown)
    dt = temps[1:] - temps[:-1]
    df = f[1:] - f[:-1]
    D = dt**2 / df  
    eta = 1./dt
    j = D*eta*df/dt
    etaopt = 1./np.sqrt(D)
    Teta = temps[:-1] + dt/2
    if True:
        print "nup", nup
        print "ndown", ndown
        print nup / (nup + ndown)
        print "diffusion", D
        print "current j", j
        with open("out.ex1", "w") as fout:
            for i in range(nreps):
                fout.write( "%g %g %g %g\n" % ( temps[i], f[i], nup[i], ndown[i]) )
            fout.write("\n\n")
            for i in range(nreps-1):
                fout.write( "%g %g %g %g %g %g %g\n" % ( Teta[i], dt[i], df[i], D[i], eta[i], j[i], etaopt[i]) )
    if False:
        import matplotlib.pyplot as plt
        plt.plot(f, temps[:])
        plt.show()
    return dt, f, D, etaopt





def getOptimalTemperatures(temps, dT, eta, nrepsnew = None, overlapmin = 0.001, Trangenew = None):
    """
    get an optimal set of temperatures from the optimal density, eta

    must be careful with the first and last temperature because eta is somewhat ill 
    defined there being that it is taken calculated from a derivative of a discrete set of data
    """
    """note: are these really the right temperatures?"""
    if nrepsnew == None:
        nrepsnew = len(temps)
    if Trangenew == None:
        Trangenew = [temps[0], temps[-1]]
    Teta = temps[:-1] + dT/2
    eta /= np.sum( (eta[1:]+eta[:-1])/2. * (Teta[1:]-Teta[:-1]))
    if True:
        #calculate the temperatures analytically from eta
        nT = 10000
        Tmin = temps[0]
        Tmax = temps[-1]
        """note: Tmin and Tmax are outside the bounds of Teta the interpolation might give strange results in these cases"""
        Tinterp = np.array( [Tmin + float(i)/(nT-1)*(Tmax-Tmin) for i in range(nT) ] )
        etainterp = np.interp(Tinterp, Teta, eta)
        etainterp /= np.sum( (etainterp[1:]+etainterp[:-1])/2. * (Tinterp[1:] - Tinterp[:-1]))
        Topt = np.zeros(nrepsnew)
        Topt[0] = temps[0]
        Topt[-1] = temps[-1]
        k = 1
        integral = 0.
        integrand =  (etainterp[1:]+etainterp[:-1])/2. * (Tinterp[1:] - Tinterp[:-1])
        cumsum = np.cumsum(integrand)
        print "cumsum", cumsum
        for i in range(1,len(Tinterp)-1):
            if cumsum[i] >= float(k) / (nrepsnew):
                Topt[k] = Tinterp[i]
                k += 1
                #print "temp ", k, Topt[k]
                if k == nrepsnew-1:
                    break
        print Topt
    else:
        #use evenly spaced temperatures as a starting point for the quench
        dT = (temps[-1] - temps[0])/ (nrepsnew-1)
        Topt = np.array( temps[0] + dT*i for i in range(nrepsnew) )
    if True:
        #use geometrically distributed temperatures as a starting point for the quench
        Tinit = np.zeros(nrepsnew)
        CTE=(np.log(Tmax/Tmin))/(nrepsnew-1)
        CTE=np.exp(CTE)
        #for J1 in range(nrepsnew):
        Tinit[:] = Tmin*CTE**range(nrepsnew)

    means, stds = readVisitsFiles()

    """
    now try doing it more analytically but including the minimum overlap 
    """
    from get_exchanges_Tdist import TCalc, Eta, Overlap, TCalcNew
    Qclass = Overlap(means, stds**2, temps, overlapmin)
    etaclass = Eta(Teta, eta, Tmin, Tmax)
    if False:
        calc = TCalc(Tmin, Tmax, nrepsnew, etaclass, Qclass, overlapmin)
        calc.findTemps()
        Toptquenched = copy.copy(calc.temps)
    if True:
        calc = TCalcNew(Trangenew[0], Trangenew[-1], nrepsnew, etaclass, Qclass, overlapmin)
        calc.findTemps()
        Toptquenched = copy.copy(calc.temps)


    if False:
        """
        the analytical result seems to mess up a bit far from the heat capacity peak.
        It gives replicas that are so far apart in temperature that they very rarely interact.
    
        Here we use a minimization to ensure that all replicas have a minimum overlap while trying
        to arrange them according to eta as well as possible
        """
        from topt_pot import TOptPot
        pot = TOptPot(eta, Teta, temps[0], temps[-1], nrepsnew, means, stds, temps, overlapmin=overlapmin)
        from pele.optimize.quench import lbfgs_py as quench
        #WARNING, the required tolerance seems to be completely system dependent
        ret = quench( Tinit, pot.getEnergyGradient, maxstep=0.1, tol = 0.005, iprint=1)
        Toptquenched = ret[0]
    
    
    if True:
        volumes = [ etaclass.integrate( Toptquenched[i], Toptquenched[i+1] ) for i in range(nrepsnew-1) ]
        overlaps = Qclass(Toptquenched)
        for i in range(nrepsnew-1):
                print "volumes, overlaps", Toptquenched[i], Toptquenched[i+1], volumes[i], overlaps[i]
        with open("exchanges.Topt", "w") as fout:
            for i in range(nrepsnew):
                fout.write( "%g %g %g\n" % ( Toptquenched[i], Topt[i], Tinit[i] ) )
            fout.write("\n\n")
                
    return Toptquenched


def getOverlaps(temps, means, stds, newtemps = None):
    from get_exchanges_Tdist import getOverlap
    if newtemps == None:
        newtemps = temps
    x = np.interp(newtemps, temps, means)
    s = np.interp(newtemps, temps, stds)
    overlaps = []
    for i in range(len(newtemps) - 1):
        q = getOverlap( x[i+1], s[i+1]**2, x[i], s[i]**2 )
        overlaps.append( q )
        print "overlap ", newtemps[i], newtemps[i+1], q




def readVisitsFiles():
    #load a Visits.his file to get the energy histograms
    dirlist = plotVisits.getDirList()
    vis = plotVisits.listVisitsFiles(dirlist[0])
    vis = vis[-1]
    fname = vis.fname
    allvisits = []
    for d in dirlist:
        visits = utils.VisitsHis( d+"/" + fname )
        allvisits.append(visits)
    means = []
    stds = []
    for v in allvisits:
        v.read()
        energies, visits = v.getHis()
        m, std = utils.wstd(energies, visits)
        means.append(m)
        stds.append(std)
    return np.array(means), np.array(stds)

def estimateNewExchangeProb(temps, newtemps):
    means, stds  = readVisitsFiles()
    if True:
        for i in range(len(temps)):
            print temps[i], means[i], stds[i] 
    print ""
    print "overlaps old temps"
    getOverlaps(temps, means, stds)
    #print ""
    #print "overlaps new temps"
    #overlaps = getOverlaps(temps, means, stds, newtemps)
    #return means, stds, overlaps






def getNupNdown(exchanges, nreps):
    """
    http://dx.doi.org/10.1088/1742-5468/2006/03/P03018

    nup =   the number of times a temperature was occupied by a replica who last visited Tmax
    ndown = the number of times a temperature was occupied by a replica who last visited Tmin
    """
    up = 1
    down = -1
    neither = 0
    state = np.ones(nreps) * neither
    state[0] = down
    state[-1] = up
    nup = np.zeros(nreps)
    ndown = np.zeros(nreps)
    timeprev = exchanges[0][0]  #this ignores the first exchange time
    for ex in exchanges:
        time = ex[0]
        rep1 = ex[1]
        rep2 = ex[2]

        #update nup and ndown
        dt = time - timeprev
        timeprev = time
        ind = np.where(state == up)[0]
        nup[ind] += dt
        ind = np.where(state == down)[0]
        ndown[ind] += dt


        #update state
        t = state[rep1]
        state[rep1] = state[rep2]
        state[rep2] = t

        #if either is at one of the boundaries, change the state
        if rep1 == 0:
            state[rep1] = down
        if rep2 == nreps-1:
            state[rep2] = up
    return nup, ndown



def main():
    ############################################################
    # read command line parameters
    ############################################################
    import getopt, sys
    import copy
    if len(sys.argv) <= 1:
        dirlist = '.'
    else:
        dirlist = copy.copy(sys.argv[1:])



    fpkl = "exchanges.pickle"
    nreps = utils.getNReplicas()
    temps = np.loadtxt("temperatures")
    print "using nreps = ", nreps
    try:
        with open(fpkl, "r") as fin:
            print "loading exchanges from pickle file", fpkl
            exall, diffusiondata = pickle.load(fin)
    except:
        print "reading exchanges from directories", dirlist
        #dirlist = [ './', '../cavity200-5-restart2' ]
        readclass = ReadExchanges( nreps, dirlist )
        exall = readclass.readAll()
        print exall[0][:3]
        #exall = read_exchanges(nreps)
        diffusiondata = getDiffusion(exall, nreps, temps)
        with open(fpkl, "w") as fout:
            pickle.dump((exall, diffusiondata), fout)

    nexch = [0 for i in range(nreps-1)]
    nexch = np.zeros( nreps-1, np.integer )
    for ex in exall:
        nexch[ int(ex[1])-1 ] += 1
    with open("exchanges.dat", "w") as fout:
        fout.write("#i i+1 nexchanges\n")
        for i in range(nreps-1):
            print "number of exchanges between %3d and %3d = %10d" % (i+1, i+2, nexch[i])
            fout.write( "%3d %3d %10d\n" % (i+1, i+2, nexch[i]) )

    if True:
        try:
            dt, f, D, etaopt = diffusiondata
            Topt = getOptimalTemperatures(temps, dt, etaopt, nrepsnew = 40, overlapmin = 0.0010, Trangenew = [temps[0], temps[-1]])
            estimateNewExchangeProb(temps, Topt)
        except ImportError:
            print "WARNING: excption raised during while analyzing replica diffusion.  is scipy not installed?"
        
    if True:
        try:
            dt, f, D, etaopt = diffusiondata
            fname = "exchanges.eta.pdf"
            print "saving figure to", fname
            Thalf = (temps[1:] + temps[:-1])/2
            pp=PdfPages(fname)
            plt.clf()
            df = f[1:] - f[:-1]
            myeta = MyEta(temps, f)
            newT = np.linspace(temps[0], temps[-1], 100)
            neweta = myeta.eta(newT)
            newdf = myeta.dfdT(newT)
            if True:
                plt.plot(temps, f, label='f')
                plt.plot(newT, myeta.finterp(newT), label='f')
                plt.plot(Thalf, etaopt / max(etaopt), label='eta optimal')
                plt.plot(newT, neweta / max(neweta), label='eta optimal')
                plt.plot(Thalf, df / max(df), label='df')
                plt.plot(newT, newdf / max(newdf), label='df')
                plt.xlabel('T')
                plt.legend()
                pp.savefig()
                plt.clf()
            if True:
                plt.plot(f, label='f')
                plt.plot(etaopt / max(etaopt), label='eta optimal')
                plt.plot(df / max(df), label='df')
                plt.xlabel('replica number')
                plt.legend()
                pp.savefig()
                plt.clf()
            pp.close()
        except ImportError:
            pass


    if True:
        pathlist = []
        roundtrips = []
        for rep in range(nreps):
            path, trips = find_path(rep, exall, nreps)
            pathlist.append( path )
            roundtrips.append( trips )
        #for path in pathlist:
            #print "path of replica", path[0]
            #for rep in path:
                #print rep
        
        for rep, trips in enumerate(roundtrips):
            print "complete trips up, down for replica %3d: %4d %4d" % ( rep, trips[0], trips[1])

        fname="exchanges.pdf"
        pp=PdfPages(fname)
        plt.plot( nexch, '-.' )
        plt.xlabel( "i" )
        plt.ylabel( "number of exchanges i -> i+1" )
        pp.savefig()
        plt.clf()
        timemax = max( [ path[-1,0] for path in pathlist ] )
        timemin = min( [ path[0,0] for path in pathlist ] )
        print "timemax, timemin", timemax, timemin
        count = 0
        for path in pathlist:
            count += 1
            path = np.array(path)
            #print path[:,1].min()
            temppath = np.array( [ temps[int(rep)] for rep in path[:,1]])
            for T in temps:
                plt.plot( [timemin, timemax], [T,T], '-k')
            plt.plot( path[:,0], temppath, '-b' )
            plt.ylim((temps[0],temps[-1]))
            plt.xlim(xmin = timemin, xmax = timemax)
            plt.ylabel("temperature slot")
            plt.xlabel("MC steps")
            plt.title("path of replica "+str(count) )
            #plt.show()
            pp.savefig()
            plt.clf()
        pp.close()

if __name__ == "__main__":
    main()
