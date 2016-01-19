import numpy as np
from numpy import log, exp
from scipy.misc import logsumexp
from scipy.optimize import fmin_l_bfgs_b
try:
    from scipy.optimize import Result
except ImportError:
    # they changed the name from Result to OptimizeResult at some point
    from scipy.optimize import OptimizeResult as Result

def dos_from_offsets(visits, log_dos_all, offsets, nodata_value=0.):
    log_dos_all = log_dos_all + offsets[:,np.newaxis]
    
    ldos = np.sum(log_dos_all * visits, axis=0)
    norm = visits.sum(0)
    ldos = np.where(norm > 0, ldos / norm, nodata_value)
    return ldos

def estimate_dos(visits, reduced_energy):
    """estimate the density of states from the histograms of bins
    
    Notes
    -----
    The density of states is proportional to visits * exp(E/T).  Therefore
    the log density of states for each replica is related to one another by an
    additive constant.  This function will find those offsets and produce a
    guess for the total density of states.
    
    """
    SMALL = 0.
    if len(visits.shape) != 2:
        raise ValueError("visits must have two dimensions")
    if visits.shape != reduced_energy.shape:
        raise ValueError("visits has a different shape from reduced_energy")
    log_dos = np.where(visits==0, SMALL, 
                       np.log(visits) + reduced_energy) 
    
    offsets = [0.]
    nreps = visits.shape[0]
    for i in xrange(1,nreps):
        # find the difference in the log density of states
        ldos_diff = log_dos[i-1,:] - log_dos[i,:]
        # weight the difference by the minimum visits in each bin
        weights = np.where(visits[i-1,:] < visits[i,:], visits[i-1,:], visits[i,:])
        new_offset = np.average(ldos_diff, weights=weights)
        offsets.append( offsets[-1] + new_offset)
    offsets = np.array(offsets)
    ldos = dos_from_offsets(visits, log_dos, offsets)

    return offsets, ldos


def calc_Cv(Tlist, binenergy, logn_E, ndof, have_data=None, k_B=1.):
    """compute the heat capacity and other data from the density of states
    
    Parameters
    ----------
    Tlist : list
        list of temperatures at which to do the computations
    binenergy : array
        The lower edge of the bins for the density of states data
    logn_E : numpy array
        the binned density of states.
    ndof : int
        number of degrees of freedom
    have_data : numpy array
        array of indices where data is available.  If have_data is None 
        then all data is assumed good
    k_B : float
        Boltzmann's constant.  Conversion factor between temperature and energy
    
    """
    nz = have_data
    
    assert(not np.any(np.isnan(logn_E[nz])))
    assert(not np.any(np.isnan(binenergy[nz])))
    
    Ebin_min = binenergy[nz].min() - 1.
    nebins = len(binenergy)

    # now calculate partition functions, energy expectation values, and Cv
    dataout = np.zeros([len(Tlist), 6])
    if abs((binenergy[-1] - binenergy[-2]) - (binenergy[-2] - binenergy[-3]) ) > 1e-7:
        print "calc_Cv: WARNING: dE is not treated correctly for uneven energy bins"

    for count, T in enumerate(Tlist):
        kBT = k_B * T
        # find expoffset so the exponentials don't blow up
        dummy = logn_E[nz] - binenergy[nz] / kBT

        lZ0 = logsumexp(dummy)
        lZ1 = logsumexp(dummy + log(binenergy[nz] - Ebin_min))
        lZ2 = logsumexp(dummy + 2. * log(binenergy[nz] - Ebin_min))

        i = nebins-1
        #if abs(binenergy[-1] - binenergy[-2]) > 1e-7:
        #    print "calc_Cv: WARNING: dE is not treated correctly for exponential bins"
        if i == nebins-1:
            dE = binenergy[i] - binenergy[i-1]
        else:
            dE = binenergy[i+1] - binenergy[i]
            
        if dE/kBT < 1e-7:
            OneMExp = -dE/kBT
        else:
            OneMExp = 1.0 - np.exp(dE/kBT)

        Eavg = ndof * kBT / 2.0 + kBT + dE / OneMExp + exp(lZ1-lZ0) + Ebin_min
        
        # correction to the Cv for finite size of bins
        bin_correction = 1.0 - dE**2 * exp(dE / kBT) / (OneMExp**2 * kBT**2)
        Cv = k_B * (ndof / 2. + (exp(lZ2-lZ0) - exp(lZ1 - lZ0)**2) / kBT**2
                    + bin_correction
                    )
        
        dataout[count,0] = T
        dataout[count,1] = lZ0
        dataout[count,2] = lZ1
        dataout[count,3] = lZ2
        dataout[count,4] = Eavg
        dataout[count,5] = Cv
        
    return dataout


def lbfgs_scipy(coords, pot, iprint=-1, tol=1e-3, nsteps=10000):
    """
    a wrapper function for lbfgs routine in scipy
    
    .. warn::
        the scipy version of lbfgs uses linesearch based only on energy
        which can make the minimization stop early.  When the step size
        is so small that the energy doesn't change to within machine precision (times the
        parameter `factr`) the routine declares success and stops.  This sounds fine, but
        if the gradient is analytical the gradient can still be not converged.  This is
        because in the vicinity of the minimum the gradient changes much more rapidly then
        the energy.  Thus we want to make factr as small as possible.  Unfortunately,
        if we make it too small the routine realizes that the linesearch routine
        isn't working and declares failure and exits.
        
        So long story short, if your tolerance is very small (< 1e-6) this routine
        will probably stop before truly reaching that tolerance.  If you reduce `factr` 
        too much to mitigate this lbfgs will stop anyway, but declare failure misleadingly.  
    """
    assert hasattr(pot, "getEnergyGradient")
    res = Result()
    res.coords, res.energy, dictionary = fmin_l_bfgs_b(pot.getEnergyGradient, 
            coords, iprint=iprint, pgtol=tol, maxfun=nsteps, factr=10.)
    res.grad = dictionary["grad"]
    res.nfev = dictionary["funcalls"]
    warnflag = dictionary['warnflag']
    #res.nsteps = dictionary['nit'] #  new in scipy version 0.12
    res.nsteps = res.nfev
    res.message = dictionary['task']
    res.success = True
    if warnflag > 0:
        print "warning: problem with quench: ",
        res.success = False
        if warnflag == 1:
            res.message = "too many function evaluations"
        else:
            res.message = str(dictionary['task'])
        print res.message
        print "    the energy is", res.energy, "the rms gradient is", np.linalg.norm(res.grad) / np.sqrt(res.grad.size), "nfev", res.nfev
        print "    X: ", res.coords
    #note: if the linesearch fails the lbfgs may fail without setting warnflag.  Check
    #tolerance exactly
    if False:
        if res.success:
            maxV = np.max( np.abs(res.grad) )
            if maxV > tol:
                print "warning: gradient seems too large", maxV, "tol =", tol, ". This is a known, but not understood issue of scipy_lbfgs"
                print res.message
    res.rms = res.grad.std()
    return res

