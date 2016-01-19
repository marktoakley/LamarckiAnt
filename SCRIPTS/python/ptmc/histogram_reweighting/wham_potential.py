import numpy as np

from pele.potentials import BasePotential

class WhamPotential(BasePotential):
    """
    the idea behind this minimization procedure is as follows 
    
    from a simulation at temperature T you find the probability of finding energy
    E is P(E,T).  We know this can be compared to the density of states n(E) as
    
        P(E,T) = n(E) exp(-E/T_i) / w_i
    
    Where w_i is a constant that is not known.  The density of 
    states is independent of temperature, so we can use it to find
    P(E) at any other temperature, or Z(T), etc.  But our estimate of n(E) from
    one temperature is not very good.  So we combine P(E,T) multiple simulations
    at different temperatures to get a better estimate of n(E).  
    Define R the log deviation for each bin from the estimate of the density of states
    
        R(E,T_i) = log(n_F(E)) - log(w_i) - log( P(E,T_i) * exp(E/T_i) )
    
    we want to make each R(E,T_i) as small as possible.  Define an "energy" function
    
        CHI2 = sum_E sum_i P(E,T_i) * |R(E,T_i)|^2
    
    Where each R(E,T_i) contributes weight proportional to P(E,T_i) to the sum to
    make sure those with better statistics are more heavily weighted.  To solve
    the problem we find the set of {n_F(E), w_i} which minimize CHI2
    """
    def __init__(self, P, reduced_energy):
        """
        To make it fit within existing minimization schemes, we need to view it as a linear problem

        nrep:  the number of replica variables, i.e. len(w_i)

        nbins: the number of bins in the histogram, e.g. len(n_F(E))

        P: = P(E,T_i) a.k.a. log(visits) a 2d array of shape( nreps, nbins).

        reduced_energy:  E/T_i  a 2d array of shape( nreps, nbins) giving the
            reduced energy of each bin

        note: this works perfectly well for 2d histograms as well.  In this case the 2d 
            histograms should be linearized
        
        """
        self.nreps, self.nbins = P.shape
        assert P.shape == reduced_energy.shape
        self.P = P
        
        if np.any(self.P < 0):
            raise ValueError("P has negative values")
        
        SMALL = 0. # this number is irrelevant, as long as it's not NaN
        self.log_n_rE = np.where(self.P==0, SMALL, 
                                 np.log(self.P) + reduced_energy)
        
    def getEnergy(self, X):
        """
        X: is the array of unknowns of length nrep + nbins
            X[0:nrep] = {w_i}         : the replica unknowns
            X[nrep:] = {log(n_F(E))}  : the bin unknowns

        R(E,T_i) = log(n_F(E)) - log(w_i) - log( P(E,T_i)*exp(E/T_i) )
    
        energy = sum_E sum_i P(E,T_i)*|R(E,T_i)|^2
        """
        wi = X[:self.nreps]
        lognF = X[self.nreps:]
        energy = np.sum( self.P * (lognF[np.newaxis,:] - wi[:,np.newaxis] - self.log_n_rE)**2 )
        return energy

    def getEnergyGradient(self, X):
        """
        X: is the array of unknowns of length nrep + nbins
            X[0:nrep] = {w_i}         : the replica unknowns
            X[nrep:] = {log(n_F(E))}  : the bin unknowns

        R(E,T_i) = log(n_F(E)) - log(w_i) - log( P(E,T_i)*exp(E/T_i) )
    
        energy = sum_E sum_i P(E,T_i)*|R(E,T_i)|^2
        """
        wi = X[:self.nreps]
        lognF = X[self.nreps:]
        R = lognF[np.newaxis,:] - wi[:,np.newaxis] - self.log_n_rE
        
        energy = np.sum( self.P * (R)**2 )
        
        gradient = np.zeros(len(X))
        gradient[:self.nreps] = -2. * np.sum( self.P * R, axis=1 )
        gradient[self.nreps:] = 2. * np.sum( self.P * R, axis=0 )
        #print np.shape(gradient)
        #print gradient
    
        return energy, gradient

    
