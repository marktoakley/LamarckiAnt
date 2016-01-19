import numpy as np
import scipy.interpolate
import copy

def getOverlap(x1, var1, x2, var2):
    """
    calculate the overlap of two gaussians
    """
    return np.exp( -(x1-x2)**2 / (var1+var2) ) / np.sqrt(np.pi) / np.sqrt(var1+var2)

class Overlap(object):
    def __init__(self, means, vars, Tmeans, overlapmin = 0.001):
        #set up interpolation of means and standard deviations for overlap calculation
        self.meaninterp = scipy.interpolate.UnivariateSpline( Tmeans, means, s=0)
        self.varinterp = scipy.interpolate.UnivariateSpline( Tmeans, vars, s=0)
    def getOverlaps(self, temps):
        x = self.meaninterp(temps)
        var = self.varinterp(temps)
        overlaps = np.zeros(len(temps)-1)
        for i in range(len(temps)-1):
            overlaps[i] = getOverlap( x[i], var[i], x[i+1], var[i+1] )
        return overlaps
    def __call__(self, temps):
        return self.getOverlaps(temps)

class Eta(object):
    def __init__(self, Teta, eta, Tmin, Tmax):
        #if eta doesn't cover the whole range [Tmin,Tmax], extend it linearly
        if Teta[0] > Tmin:
            etaTmin = eta[0] + (eta[1] - eta[0])/(Teta[1] - Teta[0]) * (Tmin - Teta[0])
            newT = np.zeros( 1+len(Teta))
            newT[0] = Tmin
            newT[1:] = Teta[:]
            neweta = np.zeros( 1+len(Teta))
            neweta[0] = etaTmin
            neweta[1:] = eta[:]
            Teta = newT
            eta = neweta
        if Teta[-1] < Tmax:
            etaTmax = eta[-1] + (eta[-1] - eta[-2])/(Teta[-1] - Teta[-2]) * (Tmax  - Teta[-1])
            newT = np.zeros( 1+len(Teta))
            newT[-1] = Tmax
            newT[:-1] = Teta[:]
            neweta = np.zeros( 1+len(Teta))
            neweta[-1] = etaTmax
            neweta[:-1] = eta[:]
            Teta = newT
            eta = neweta
        self.eta = eta
        self.Teta = Teta
        print "after interpolation Teta", Teta
        print "after interpolation Teta", eta
        #import matplotlib.pyplot as plt
        #plt.plot( self.Teta, self.eta, '-x')
        #plt.show()
        #exit(1)
        self.etainterp = scipy.interpolate.UnivariateSpline(self.Teta, self.eta, s=0)
        self.etanorm = self.etainterp.integral(Tmin, Tmax)
    def getEta(self, T):
        """
        interpolate to find the density at a given temperature
        """
        #return np.interp(T, self.Teta, self.eta)
        return self.etainterp(T) / self.etanorm
    def __call__(self, T):
        return self.getEta(T)
    def integrate(self, T1, T2):
        return self.etainterp.integral(T1, T2) / self.etanorm


class TCalc(object):
    """
    return temperatures distributed as evenly as possible on the interval
    Tmin, Tmax according to the distribution eta, subject to the constraint
    overlap(T[i], T[i+]) > overlapmin
    """
    def __init__(self, Tmin, Tmax, nreps, eta, getOverlaps, overlapmin):

        self.Tmin = Tmin
        self.Tmax = Tmax
        self.nreps = nreps
        self.temps = np.ones(self.nreps)
        self.temps[0] = Tmin
        self.temps[-1] = Tmax
        self.dtemps = np.zeros(self.nreps-1)

        self.overlapmin = overlapmin

        self.count = 0
        self.eta = eta
        self.getOverlaps = getOverlaps

        self.fixed = [False for i in range(self.nreps-1) ] #contains the temperature intervals that are subject to the overlap constraint 
        self.fixed = np.array(self.fixed)
    
    def getVolume(self):
        """
        return volume per interval
        """
        self.setTemps()
        if not any(self.fixed):
            return 1./(self.nreps-1)
        V = 0.
        for i in range(self.nreps-1):
            if not self.fixed[i]: continue
            vol = self.eta.integrate(self.temps[i], self.temps[i+1])
            #print "integral", self.temps[i], self.temps[i+1], vol
            V += vol
        nfixed = sum(self.fixed)
        print "nfixed", nfixed, V
        return (1.-V)/ (self.nreps-1 - nfixed)
            
    
    def placeTemps(self):
        """
        starting from Tmin, place temps such that the integral over eta of the temperature
        intervals is equal.  
        """
        fixed  = self.fixed
        temps = self.temps
        vol = self.getVolume()
        dT = 1e-3
        T = self.Tmin
        irep = 1
        while irep < self.nreps-1:
            T += dT
            if T >= self.Tmax:
                #raise Exception("T>Tmax. irep=%d" %(irep) )
                print "warning: T>Tmax. irep=%d" %(irep)
                fixed[irep:] = False
                if False:
                    #randomly unfix an interval
                    nfixed = sum(fixed)
                    irand = np.random.random_integers(1, nfixed)
                    count = 0
                    for i in range(len(fixed)):
                        if fixed[i]:
                            count+=1
                            if count == irand:
                                print "randomly unfixing interval", i
                                fixed[i] = False
                                break
                if False:
                    #unfix the interval with the largest integrated volume
                    vmax = 0
                    imax = 0
                    for i in range(len(fixed)):
                        if fixed[i]:
                            vol = self.eta.integrate(self.temps[i], self.temps[i+1])
                            if vol > vmax:
                                vmax = vol
                                imax = i
                    print "unfixing interval", imax, "with volume", vmax
                    fixed[imax] = False
                return False
            
            v = self.eta.integrate(self.temps[irep-1], T)
            if v >= vol:
                #place temperature
                self.temps[irep] = T
                if fixed[irep-1]:
                    #an interval has become unfixed.  mark it as unfixed and abort
                    print "unfixing interval", irep-1
                    fixed[irep-1] = False
                    fixed[irep:] = False #unfix all the ones we didn't get to
                    return False

                irep += 1
                continue
            
            q = self.getOverlaps( np.array([self.temps[irep-1], T]) )
            if q < self.overlapmin:
                self.temps[irep] = T
                if not fixed[irep-1]:
                    #we have a newly fixed interval.  mark it as fixed and abort 
                    print "fixing interval", irep-1
                    fixed[irep-1] = True
                    #fixed[irep:] = False #unfix all the ones we didn't get to
                    return False
                irep += 1
                continue
        
        """
        we've made it to the end, check if the overlap with Tmax is ok
        """
        q = self.getOverlaps( np.array([self.temps[-2], self.temps[-1]]) )
        if q < self.overlapmin:
            fixed[-1] = True
            return False

        return True
            
    def findTemps1(self):
        
        while True:
            ret = self.placeTemps()
            print "T = ", self.temps
            if ret:
                break

    def findTemps(self):        
        while True:
            ret = self.placeTempsInternal()
            self.setTemps()
            print "T = ", self.temps
            if ret:
                break


    def dTtoTemp(self, i):
        if i == self.nreps-1: return self.Tmax
        return self.Tmin + np.sum(self.dtemps[:i])
    def setTemps(self):
        for i in range(self.nreps):
            self.temps[i] = self.dTtoTemp(i)


    def placeTempsInternal(self):
        """
        starting from Tmin, place temps such that the integral over eta of the temperature
        intervals is equal.  
        """
        #self.temps = np.array([ self.dTtoTemp(i) for i in range(self.ntemps) ])
        self.setTemps()
        fixed  = self.fixed
        temps = self.temps
        vol = self.getVolume()
        dT = 1e-3
        T = self.Tmin
        irep = 1
        while irep < self.nreps-1:
            T += dT
            if T >= self.Tmax:
                #raise Exception("T>Tmax. irep=%d" %(irep) )
                print "warning: T>Tmax. irep=%d" %(irep)
                fixed[irep:] = False
                if False:
                    #randomly unfix an interval
                    nfixed = sum(fixed)
                    irand = np.random.random_integers(1, nfixed)
                    count = 0
                    for i in range(len(fixed)):
                        if fixed[i]:
                            count+=1
                            if count == irand:
                                print "randomly unfixing interval", i
                                fixed[i] = False
                                break
                if False:
                    #unfix the interval with the largest integrated volume
                    vmax = 0
                    imax = 0
                    for i in range(len(fixed)):
                        if fixed[i]:
                            vol = self.eta.integrate(self.dTtoTemp(i), self.dTtoTemp(i+1))
                            if vol > vmax:
                                vmax = vol
                                imax = i
                    print "unfixing interval", imax, "with volume", vmax
                    fixed[imax] = False
                return False
            
            v = self.eta.integrate(self.dTtoTemp(irep- 1), T)
            if v >= vol:
                #place temperature
                self.dtemps[irep-1] = T - self.dTtoTemp(irep-1)
                if fixed[irep-1]:
                    #an interval has become unfixed.  mark it as unfixed and abort
                    print "unfixing interval", irep-1
                    fixed[irep-1] = False
                    fixed[irep-1:] = False #unfix all the ones we didn't get to
                    return False

                irep += 1
                continue
            
            q = self.getOverlaps( np.array([self.dTtoTemp(irep-1), T]) )
            if q < self.overlapmin:
                self.dtemps[irep-1] = T - self.dTtoTemp(irep-1)
                if not fixed[irep-1]:
                    #we have a newly fixed interval.  mark it as fixed and abort 
                    print "fixing interval", irep-1, self.dtemps[irep-1], T, self.dTtoTemp(irep-1)
                    fixed[irep-1] = True
                    #fixed[irep:] = False #unfix all the ones we didn't get to
                    return False
                irep += 1
                continue
        
        """
        we've made it to the end, check if the overlap with Tmax is ok
        """
        #q = self.getOverlaps( np.array([self.temps[-2], self.temps[-1]]) )
        q = self.getOverlaps( np.array( [self.dTtoTemp(self.nreps-1), self.dTtoTemp(self.nreps)] ) )
        if q < self.overlapmin:
            fixed[-1] = True
            return False
        #if the volume of the last interval is not good then run it one more time
        v = self.eta.integrate( self.dTtoTemp(self.nreps-1), self.dTtoTemp(self.nreps) )
        if np.abs((v - vol)/vol) > 1e-3:
            return False
            


        return True
  
from pele.potentials.potential import potential as basepot
class TCalcNew(basepot):
    """
    return temperatures distributed as evenly as possible on the interval
    Tmin, Tmax according to the distribution eta, subject to the constraint
    overlap(T[i], T[i+]) > overlapmin
    """
    def __init__(self, Tmin, Tmax, nreps, eta, getOverlaps, overlapmin):

        self.Tmin = Tmin
        self.Tmax = Tmax
        self.nreps = nreps
        self.temps = np.ones(self.nreps)
        self.temps[0] = Tmin
        self.temps[-1] = Tmax
        self.dtemps = np.zeros(self.nreps-1)

        self.overlapmin = overlapmin

        self.count = 0
        self.eta = eta
        self.getOverlaps = getOverlaps

        self.fixed = [False for i in range(self.nreps-1) ] #contains the temperature intervals that are subject to the overlap constraint 
        self.fixed = np.array(self.fixed)
    
            
    
    def placeTemps(self, vol):
        """
        starting from Tmin, place temps such that the integral over eta of the temperature
        intervals is equal.  
        """
        fixed  = self.fixed
        temps = self.temps
        #vol = self.getVolume()
        dT = 1e-4
        T = self.Tmin
        irep = 1
        while irep < self.nreps-1:
            T += dT
            if T >= self.Tmax:
                break
            
            v = self.eta.integrate(self.temps[irep-1], T)
            if v >= vol:
                #place temperature
                self.temps[irep] = T
                irep += 1
                continue
            
            q = self.getOverlaps( np.array([self.temps[irep-1], T]) )
            if q < self.overlapmin:
                self.temps[irep] = T
                irep += 1
                continue
        

        return irep
    
    
    
    def increaseVolume(self, volume):
        if self.status == self.decreased:
            if self.reversed == self.max_reversed:
                self.dV /= 10
            self.reversed += 1
        self.status = self.increased
        return volume + self.dV
    
    def decreaseVolume(self, volume):
        if self.status == self.increased:
            if self.reversed == self.max_reversed:
                self.dV /= 10
            self.reversed += 1
        self.status = self.decreased
        return volume - self.dV

        
        
         
    def findTemps(self):
        self.increased = 1
        self.decreased = -1
        self.status = 0
        self.reversed = 0
        self.max_reversed = 10
        volume = 1./(self.nreps-1)
        dV = volume / 100
        self.dV = dV
        for count in range(200):
            nplaced = self.placeTemps(volume)
            if nplaced < self.nreps-1:
                #not enough placed.  volume is too big
                volume = self.decreaseVolume(volume)
            elif nplaced == self.nreps-1:
                q = self.getOverlaps( np.array([self.temps[-2], self.Tmax]) )
                v = self.eta.integrate(self.temps[-2], self.Tmax)
                #the last rep was placed. check the volume of the last interval
                if q < self.overlapmin:
                    volume = self.increaseVolume(volume)
                elif v > volume:
                    #increase volume a bit
                    volume = self.increaseVolume(volume)
                else:
                    volume = self.decreaseVolume(volume)
                print count, "volume", volume, nplaced, q, v
            print "T = ", self.temps
            if self.reversed > self.max_reversed*2:
                break

