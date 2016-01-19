import unittest
from itertools import izip
import numpy as np

import matplotlib.pyplot as plt

from histogram_reweighting1d import Wham1d
from histogram_reweighting import wham_utils

def get_temps(Tmin, Tmax, nreplicas):
    """
    set up the temperatures
    distribute them exponentially
    """
    #dT = (Tmax - Tmin) / (nreplicas-1)
    CTE = np.exp( np.log( Tmax / Tmin ) / (nreplicas-1) )
    Tlist = [Tmin* CTE**i for i in range(nreplicas)]
    return Tlist


class HarmonicOscillator(object):
    def __init__(self, d):
        self.d = d
    
    def random_energy(self, T):
        x = np.random.normal(size=self.d)
        x2 = x.dot(x)
        E = T * x2 / 2
        return E

    def random_energies(self, N, T):
        x = np.random.normal(size=self.d * N)
        x = x.reshape([N,self.d])
        x2 = np.sum(x**2, axis=1)
        assert x2.size == N
        Elist = x2 * T / 2
        return Elist

    def make_analytic_histogram(self, N, binenergy, T):
        assert self.d > 2
        binenergy = np.append(binenergy, binenergy[-1] + (binenergy[-1] - binenergy[-2]))
        ecenters = (binenergy[1:] + binenergy[:-1]) / 2
        assert np.all(ecenters > 0)
        logP = float(self.d-2) / 2 * np.log(ecenters) - ecenters / T
        
             
        
        # now normalize it to N
        import scipy
        logsumP = scipy.misc.logsumexp(logP)
        logP = logP + np.log(N) - logsumP
        visits = np.exp(logP)
#        visits = np.round(visits)
        if True:
            e1 = np.average(ecenters, weights=visits) 
            e2 = np.average(ecenters**2, weights=visits)
            print "Cv from analytic histogram", (e2 - e1**2) / T**2 + float(self.d)/2
        
#        if True:
#            plt.clf()
#            plt.plot(binenergy[:-1], visits.transpose())
#            plt.show()
        
        return visits 
        

    def random_histogram(self, N, binenergy, T):
        Elist = self.random_energies(N, T)
    
        evar = np.var(Elist)
        print "Cv computed from energy variance", evar / T**2 + float(self.d)/2
        
        binenergy = np.append(binenergy, binenergy[-1]+(binenergy[-1] - binenergy[-2]))
        
        counts, bins = np.histogram(Elist, binenergy)
        
        if True:
            e = binenergy[:-1]
#            plt.plot(e, counts)
            plt.plot(e, np.log(counts) + e / T )
        return counts
    
    def make_visits(self, Tlist, binenergy, N, random=True):
#        binenergy = np.linspace(0, 40, 1000)
#        visits = np.zeros([len(Tlist), len(binenergy)])
        visits = []
        for i, T in enumerate(Tlist):
            if random:
                counts = self.random_histogram(N, binenergy, T)
            else:
                counts = self.make_analytic_histogram(N, binenergy, T)
            visits.append(counts)
#        plt.show()

        if False:
            plt.clf()
            plt.title("before numpy array")
            for T, counts in izip(Tlist, visits):
                log_nE = np.log(counts) + binenergy / T
                plt.plot(binenergy, log_nE)
            plt.show()

        visits = np.array(visits)


        if False:
            plt.clf()
            print visits.shape, binenergy.shape, Tlist.shape
            log_nET = np.log(visits) + binenergy[np.newaxis, :] / Tlist[:,np.newaxis] #+ wham.w_i_final[:,np.newaxis]
            plt.plot(binenergy, np.transpose(log_nET))
            plt.show()
#        raise Exception

        
        return visits
    
    
    def make_analytic_dos(self, binenergy):
        """make a binned density of states
        
            dos = E**(d/2)
        """
        assert self.d > 2
        binenergy = np.append(binenergy, binenergy[-1] + (binenergy[-1] - binenergy[-2]))
        ecenters = (binenergy[1:] + binenergy[:-1]) / 2
        assert np.all(ecenters > 0)
        log_dos = float(self.d-2) / 2 * np.log(ecenters)
        return log_dos


class TestHistogramReweighting(unittest.TestCase):
    def setUp(self):
        self.d = 3
        self.N = 100000
        self.Tlist = get_temps(.2, 1.6, 8)
#        self.Tlist = [2.000000000000000111e-01,
#                2.691800385264712658e-01,
#                3.622894657055626966e-01,
#                4.876054616817901977e-01,
#                6.562682848061104357e-01,
#                8.832716109390498227e-01,
#                1.188795431309558781e+00,
#                1.600000000000000089e+00]
        self.Tlist = np.array(self.Tlist)
        self.binenergy = np.linspace(0, 20, 1000)

    def test1(self):
        np.random.seed(0)
        ho = HarmonicOscillator(self.d)
        self.visits = ho.make_visits(self.Tlist, self.binenergy, self.N, random=True)
        assert self.visits.shape == (len(self.Tlist), len(self.binenergy))
        visits = self.visits
        binenergy = self.binenergy
        Tlist = self.Tlist
        if False:
            plt.clf()
            print visits.shape, binenergy.shape, Tlist.shape
            log_nET = np.log(visits) + binenergy[np.newaxis, :] / Tlist[:,np.newaxis] #+ wham.w_i_final[:,np.newaxis]
            plt.plot(binenergy, np.transpose(log_nET))
#            for T, counts in izip(Tlist, visits):
#                log_nE = np.log(counts) + binenergy / T
#                plt.plot(binenergy, log_nE)
            plt.show()
        
        wham = Wham1d(Tlist, binenergy, visits.copy())
        wham.minimize()
        cvdata = wham.calc_Cv(3, Tlist=Tlist)
#        print cvdata.shape
#        print cvdata
#        print cvdata
#        print "Cv values", cvdata[:,5]
        
        for cv in cvdata[:,5]:
            self.assertAlmostEqual(cv, self.d, delta=.3)

        if False:
            plt.clf()
            plt.plot(Tlist, cvdata[:,5])
            plt.show()
        
        
        
        if False:
            plt.clf()
            plt.plot(binenergy, wham.logn_E)
            plt.show()
        if False:
            plt.clf()
            for r, T in enumerate(Tlist):
                v = visits[r,:]
    #            plot(bin,v)
            plt.plot(binenergy, np.log(np.transpose(visits)))
            plt.show()
        if False:
            plt.clf()
            log_nET = np.log(visits) + binenergy[np.newaxis, :]  / Tlist[:,np.newaxis] + wham.w_i_final[:,np.newaxis]
            plt.plot(binenergy, np.transpose(log_nET))
            plt.plot(binenergy, wham.logn_E)
            plt.show()
        if False:
            plt.clf()
            plt.plot(binenergy, wham.logn_E, 'k', lw=2)
            newlogn_E = wham_utils.dos_from_offsets(Tlist, binenergy, visits,
                                                    wham.w_i_final)
            plt.plot(binenergy, newlogn_E, '--r', lw=.5)
            plt.show()
    
    def test_analytic(self):
        ho = HarmonicOscillator(self.d)
        visits = ho.make_visits(self.Tlist, self.binenergy, self.N, random=False)
        
        wham = Wham1d(self.Tlist, self.binenergy, visits)
        wham.minimize()
        cvdata = wham.calc_Cv(3, Tlist=self.Tlist)
        
        for cv in cvdata[:,5]:
            self.assertAlmostEqual(cv, self.d, delta=.01)

    
    def test2(self):
        ho = HarmonicOscillator(self.d)
        self.visits = ho.make_visits(self.Tlist, self.binenergy, self.N, random=True)
        assert self.visits.shape == (len(self.Tlist), len(self.binenergy))
        self.reduced_energy = self.binenergy[np.newaxis,:] / (self.Tlist[:,np.newaxis])

        wham_utils.estimate_dos(self.visits, self.visits)
            

if __name__ == "__main__":
    unittest.main()
