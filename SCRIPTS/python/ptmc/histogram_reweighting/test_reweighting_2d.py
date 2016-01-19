import unittest
import numpy as np

from matplotlib import pyplot as plt

import test_reweighting
from histogram_reweighting.histogram_reweighting2d import Wham2d

class HarmonicOscillatorQ(test_reweighting.HarmonicOscillator):
    """add routines to build a 2d visits histogram"""
    
    
    def random_histogram_2d(self, N, binenergy, binq, T):
        binenergy = np.append(binenergy, binenergy[-1] + (binenergy[-1] - binenergy[-2]))
        binq = np.append(binq, binq[-1] + (binq[-1] - binq[-2]))

        x = np.random.normal(size=self.d * N)
        x = x.reshape([N,self.d])
        x2 = np.sum(x**2, axis=1)
        assert x2.size == N
        Elist = x2 * T / 2
        xcoord = x[:,0] # the first coordinate
        
        evar = np.var(Elist)
        print "Cv computed from energy variance", evar / T**2 + float(self.d)/2
        
        counts, bins1, bins2 = np.histogram2d(Elist, xcoord, bins=[binenergy, binq])
        
        if False:
            c = np.where(counts==0, np.nan, counts)
            plt.imshow(c)
            plt.show()
            
            v = counts.sum(1)
            plt.clf()
            plt.plot(binenergy[:-1], v)
            plt.show()

        return counts

    def make_visits_2d(self, Tlist, binenergy, binq, N, random=True):
        visits = []
        for i, T in enumerate(Tlist):
            if random:
                counts = self.random_histogram_2d(N, binenergy, binq, T)
            else:
                counts = self.make_analytic_histogram_2d(N, binenergy, binq, T)
            visits.append(counts)

        if False:
            plt.clf()
            plt.title("before numpy array")
            for T, counts in izip(Tlist, visits):
                log_nE = np.log(counts) + binenergy / T
                plt.plot(binenergy, log_nE)
            plt.show()

        visits = np.array(visits)
        print visits.shape,  (len(Tlist), len(binenergy), len(binq))
        assert visits.shape == (len(Tlist), len(binenergy), len(binq))


        if False:
            plt.clf()
            print visits.shape, binenergy.shape, Tlist.shape
            log_nET = np.log(visits) + binenergy[np.newaxis, :] / Tlist[:,np.newaxis] #+ wham.w_i_final[:,np.newaxis]
            plt.plot(binenergy, np.transpose(log_nET))
            plt.show()
#        raise Exception

        
        return visits
        
        
class TestReweighting2d(unittest.TestCase):
    def setUp(self):
        self.d = 3
        self.N = 100000
        self.Tlist = test_reweighting.get_temps(.2, 1.6, 8)
        self.Tlist = np.array(self.Tlist)
        self.binenergy = np.linspace(0, 20, 100)
        self.binq = np.linspace(0, 7, 20)

    def test(self):
        Tlist = self.Tlist
        binenergy = self.binenergy
        binq = self.binq
        ho = HarmonicOscillatorQ(self.d)
        visits2d = ho.make_visits_2d(self.Tlist, self.binenergy, self.binq, self.N, random=True)
        
        wham = Wham2d(Tlist, binenergy, binq, visits2d)
        wham.minimize()
        cvdata = wham.calc_Cv(self.d, Tlist=self.Tlist)
        print cvdata[:,5]
        for cv in cvdata[:,5]:
            self.assertAlmostEqual(cv, self.d, delta=.7)
        
#    def test_qavg(self):
#        ho = HarmonicOscillatorQ(self.d)
#        visits2d = ho.make_visits_2d(self.Tlist, self.binenergy, self.binq, self.N, random=True)
#        wham = Wham2d(self.Tlist, self.binenergy, self.binq, visits2d)
#        wham.minimize()
#        
#        qavg = wham.calc_qavg(self.Tlist)

            
        
        
if __name__ == "__main__":
    unittest.main()
