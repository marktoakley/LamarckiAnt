import numpy as np
import scipy
from numpy import exp, log

import wham_utils






class Wham2d(object):
    """ class to combine 2d histograms of energy E and order parameter q at
    multiple temperatures into one best estimate for the histogram
  
    input will be:
    filenames : list of filenames where the data can be found
    Tlist # Tlist[k] is the temperature of the simulation in filenames[k] 
  
    binenergy = zeros(nebins, float64) #lower edge of bin energy
    binq =      zeros(nqbins, float64) #lower edge of q bins
    visits2d =  zeros([nrep,nebins,nqbins], integer) #2d histogram of data
    """
    #=============================================================================================
    # Constructor.
    #=============================================================================================
    def __init__(self, Tlist, binenergy, binq, visits2d, k_B=1., verbose=True):
        if visits2d.shape != (len(Tlist), len(binenergy), len(binq)):
            raise ValueError("Visits2d has the wrong shape")
    
        #define some parameters
        self.k_B = k_B
        self.LOGMIN = -1e200
        self.verbose = verbose
    
        self.nrep = len(Tlist)
        self.nebins = len(binenergy)
        self.nqbins = len(binq)
        self.Tlist = np.array(Tlist)
        self.binenergy = np.array(binenergy)
        self.binq = np.array(binq)
        self.visits2d = np.array(visits2d, dtype = np.integer)
        
    def minimize(self):
        #shape(visits2d) is now (nqbins, nebins, nreps)
        #we need it to be (nreps, nqbins*nebins)
        #first reorder indices
        nreps = self.nrep
        nebins = self.nebins
        nqbins = self.nqbins
        nbins = self.nebins * self.nqbins
        
        reduced_energy = np.zeros([nreps, nebins, nqbins])
        for j in range(self.nqbins):
            reduced_energy[:,:,j] = self.binenergy[np.newaxis,:] / (self.Tlist[:,np.newaxis]*self.k_B)
                    
        visits = self.visits2d
        visits = np.reshape( visits, [nreps, nbins ]) 
        reduced_energy = np.reshape( reduced_energy, [nreps, nbins])           
        self.logP = np.where( visits != 0, np.log( visits.astype(float) ), 0 )

        
        from wham_potential import WhamPotential
        whampot = WhamPotential(visits, reduced_energy )
        
        nvar = nbins + nreps
        if False:
            X = np.random.rand(nvar)
        else:
            # estimate an initial guess for the offsets and density of states
            # so the minimizer converges more rapidly
            offsets_estimate, log_dos_estimate = wham_utils.estimate_dos(visits,
                                                                         reduced_energy)
            X = np.concatenate((offsets_estimate, log_dos_estimate))
            assert X.size == nvar

        E0, grad = whampot.getEnergyGradient(X)
        rms0 = np.linalg.norm(grad) / np.sqrt(grad.size)


        from wham_utils import lbfgs_scipy
        ret = lbfgs_scipy(X, whampot)
#        print "initial energy", whampot.getEnergy(X)
#        try: 
#            from pele.optimize import mylbfgs as quench
#            ret = quench(X, whampot, iprint=10, maxstep = 1e4)
#        except ImportError:
#            from pele.optimize import lbfgs_scipy as quench
#            ret = quench(X, whampot)            

        if self.verbose:
            print "chi^2 went from %g (rms %g) to %g (rms %g) in %d iterations" % (
                E0, rms0, ret.energy, ret.rms, ret.nfev)

        
        #self.logn_Eq = zeros([nebins,nqbins], float64)
        X = ret.coords
        self.logn_Eq = X[nreps:]
        self.w_i_final = X[:nreps]
        
        self.logn_Eq = np.reshape(self.logn_Eq, [nebins, nqbins])
        self.logn_Eq = np.where( self.visits2d.sum(0) == 0, self.LOGMIN, self.logn_Eq )
        
        #renormalize logn_Eq
        #self.allzero2dind = np.where(self.visits2d.sum(2) == 0)
        #self.notzero2dind = np.where(self.visits2d.sum(2) != 0)
        #print self.allzero2dind
        #self.logn_Eq -= np.min(self.logn_Eq[ self.notzero2dind[1], self.notzero2dind[0] ])

        #self.logn_Eq[self.allzero2dind] = 0



  

  
    def calc_Fq(self, TRANGE = []):
        self.allzero2dind = np.where(self.visits2d.sum(0) == 0)

  
        #put some variables in this namespace
        nebins=self.nebins
        nqbins=self.nqbins
        binenergy=self.binenergy
        visits2d=self.visits2d
        logn_Eq=self.logn_Eq
    
        if len(TRANGE) == 0:
            NTEMP = 5 # number of temperatures to calculate expectation values
            TMAX = self.Tlist[-1]
            TMIN = self.Tlist[0]
            TINT=(TMAX-TMIN)/(NTEMP-1)
            TRANGE = [ TMIN + i*TINT for i in range(NTEMP) ]
    
        #find the occupied bin with the minimum energy
        EREF=0
        for i in range(nebins):
            if visits2d[:,i,:].sum() > 0:
                EREF = binenergy[i]
                break
    
        self.nodataq = np.where((visits2d.sum(0).sum(0)) == 0)
    
        #now calculate P(q,T)
        # P(q,T) = sum_E n(E,q)*exp(-E/T)  
        #TRANGE=range(1,9)
        self.F_q = np.zeros([nqbins,len(TRANGE)])
        F_q = self.F_q
        logP_Eq = np.zeros([nebins,nqbins])
        logP_q = np.zeros(nqbins)
        for n in range(len(TRANGE)):
            T=TRANGE[n]
            for i in range(nebins):
                logP_Eq[i,:] = logn_Eq[i,:]-(binenergy[i] - EREF)/(self.k_B*T)
      
            logP_Eq[self.allzero2dind[0], self.allzero2dind[1]] = self.LOGMIN
            expoffset = np.nanmax(logP_Eq)
            #print "T expoffset ", T, expoffset
            logP_Eq -= expoffset
            #P_q = np.exp(logP_Eq).sum(0)
            # sum over the energy
            for j in range(nqbins):
                logP_q[j] = scipy.misc.logsumexp( logP_Eq[:,j] )
            logP_q[self.nodataq] = np.NaN
            F_q[:,n] = -self.k_B*T*logP_q[:]
            fmin = np.nanmin(F_q[:,n])
            F_q[:,n] -= fmin
    
        return TRANGE,F_q
  
    def calc_qavg(self, TRANGE = []):
        """calculate the average q as a function of temperature"""
        #put some variables in this namespace
        nebins=self.nebins
        nqbins=self.nqbins
        binenergy=self.binenergy
        binq=self.binq
        visits2d=self.visits2d
        logn_Eq=self.logn_Eq
    
        if len(TRANGE) == 0:
            NTEMP = 100 # number of temperatures to calculate expectation values
            TMAX = self.Tlist[-1]
            TMIN = self.Tlist[0]
            TINT=(TMAX-TMIN)/(NTEMP-1)
            TRANGE = [ TMIN + i*TINT for i in range(NTEMP) ]
    
        #find the ocupied bin with the minimum energy
        EREF=0
        for i in range(nebins):
            if visits2d[:,i,:].sum() > 0:
                EREF = binenergy[i]
                break
    
        #don't need to recalculate it
        #self.nodataq = where((visits2d.sum(2).sum(0)) == 0)
    
        #calculate the mean q at each temperature
        self.qavg = np.zeros(len(TRANGE))
    
        #now calculate P(q,T)
        # P(q,T) = sum_E n(E,q)*exp(-E/T)  
        #TRANGE=range(1,9)
        logP_Eq = np.zeros([nebins,nqbins])
        logP_q = np.zeros(nqbins)
        for n in range(len(TRANGE)):
            T=TRANGE[n]
            for i in range(nebins):
                logP_Eq[i,:] = logn_Eq[i,:]-(binenergy[i] - EREF)/(self.k_B*T)
      
            logP_Eq[self.allzero2dind[0], self.allzero2dind[1]] = self.LOGMIN
            expoffset = logP_Eq.max()
            #print "T expoffset ", T, expoffset
            logP_Eq -= expoffset
            #P_q = np.exp(logP_Eq).sum(0)
            # sum over the energy
            for j in range(nqbins):
                logP_q[j] = scipy.misc.logsumexp(logP_Eq[:,j])
            logP_q[self.nodataq] = np.NaN
      
            #find mean q
            qmin = min(binq)
            qmin -= 0.1
            lqavg = -1.0e30
            lnorm = -1.0e30
            for i in range(0,nqbins): 
                if not np.isnan(logP_q[i]):
                    lnorm = np.logaddexp(lnorm, logP_q[i]) 
                    lqavg = np.logaddexp(lqavg, logP_q[i] + log(binq[i] - qmin))
            self.qavg[n] = exp(lqavg - lnorm) + qmin
            #print lqavg
    
        return TRANGE, self.qavg
  
  
    def calc_Cv(self, ndof, Tlist=None, ntemp=100):
        visits1d = self.visits2d.sum(2)
        have_data = np.where(visits1d.sum(0) > 0)[0]
        logn_E = scipy.misc.logsumexp(self.logn_Eq, axis=1)
#        logn_E = np.zeros(nebins)
#        for i in range(nebins):
#            logn_E[i] = scipy.misc.logsumexp(self.logn_Eq[i,:])

        if Tlist is None:
            Tlist = np.linspace(self.Tlist[0], self.Tlist[-1], ntemp)
        return wham_utils.calc_Cv(Tlist, self.binenergy, logn_E, ndof,
                                  have_data=have_data, k_B=self.k_B)

