import numpy as np #to access np.exp() not built int exp

from histogram_reweighting.wham_potential import WhamPotential
from histogram_reweighting import wham_utils

class Wham1d(object):
    """Combine 1d histograms of energy at multiple temperatures into one density of states

    Parameters
    ----------
    Tlist : list of floats
        list of temperatures
    binenergy : list of floats
        the lower edges of the energy bins
    visits1d : ndarray
        two dimensional array where visits1d[r,e] is the histogram value
        for replica r at energy bin e.
    k_B : float
        the Boltzmann constant

    """
    def __init__(self, Tlist, binenergy, visits1d, k_B=1., verbose=True):
        if visits1d.shape != (len(Tlist), len(binenergy)):
            raise ValueError("visits1d has the wrong shape")
        #define some parameters
        self.k_B = 1.

        self.nrep = len(Tlist)
        self.nebins = len(binenergy)
        self.Tlist = np.asarray(Tlist)
        self.binenergy = np.asarray(binenergy)
        self.visits1d = np.asarray(visits1d)
        
        self.verbose = verbose
        

    def minimize(self):
        """compute the best estimate for the density of states"""
        nreps = self.nrep
        nbins = self.nebins
        visitsT = (self.visits1d)
        #print "min vis", np.min(visitsT)
        #print "minlogp", np.min(self.logP)
        self.reduced_energy = self.binenergy[np.newaxis,:] / (self.Tlist[:,np.newaxis] * self.k_B)
        
        self.whampot = WhamPotential(visitsT, self.reduced_energy)
        
        if False:
            X = np.random.rand( nreps + nbins )
        else:
            # estimate an initial guess for the offsets and density of states
            # so the minimizer converges more rapidly
            offsets_estimate, log_dos_estimate = wham_utils.estimate_dos(self.visits1d,
                                                                         self.reduced_energy)
            X = np.concatenate((offsets_estimate, log_dos_estimate))

        E0, grad = self.whampot.getEnergyGradient(X)
        rms0 = np.linalg.norm(grad) / np.sqrt(grad.size)
        
        try:
            from pele.optimize import lbfgs_cpp as quench
            if self.verbose:
                print "minimizing with pele lbfgs"
            ret = quench(X, self.whampot, tol=1e-3, maxstep=1e4, nsteps=10000)
        except ImportError:
            from wham_utils import lbfgs_scipy
            if self.verbose:
                print "minimizing with scipy lbfgs"
            ret = lbfgs_scipy(X, self.whampot, tol=1e-3, nsteps=10000)
        #print "quench energy", ret.energy
        
        if self.verbose:
            print "chi^2 went from %g (rms %g) to %g (rms %g) in %d iterations" % (
                E0, rms0, ret.energy, ret.rms, ret.nfev)
        
        X = ret.coords
        self.logn_E = X[nreps:]
        self.w_i_final = X[:nreps]
        

    def calc_Cv(self, ndof, Tlist=None, ntemp=100):
        if Tlist is None:
            Tlist = np.linspace(self.Tlist[0], self.Tlist[-1], ntemp)
        have_data = np.where(self.visits1d.sum(0) > 0)[0]
        return wham_utils.calc_Cv(Tlist, self.binenergy, self.logn_E, ndof,
                                  have_data=have_data, k_B=self.k_B)
                                  

