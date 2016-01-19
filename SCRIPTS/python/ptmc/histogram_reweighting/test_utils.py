import unittest

import numpy as np

from histogram_reweighting import wham_utils
import test_reweighting


class TestCvCalc(unittest.TestCase):
    def test(self):
        d = 3
        ho = test_reweighting.HarmonicOscillator(d)
        binenergy = np.linspace(0, 40, 200)
        Tlist = np.linspace(.2, 1.6, 8)
        log_dos = ho.make_analytic_dos(binenergy)
        
        cvdata = wham_utils.calc_Cv(Tlist, binenergy, log_dos, d)
        cv = cvdata[:,5]
        self.assertLess(np.max(np.abs(cv - d)), 0.05)
        
    def test_with_zeros(self):
        d = 3
        ho = test_reweighting.HarmonicOscillator(d)
        binenergy = np.linspace(0, 40, 200)
        Tlist = np.linspace(.2, 1.6, 8)
        log_dos = ho.make_analytic_dos(binenergy)
        have_data = np.array(range(len(binenergy)-100))
#        log_dos[0] = 0.
        log_dos[-100:] = 0.
        
        cvdata = wham_utils.calc_Cv(Tlist, binenergy, log_dos, d, have_data=have_data)
        cv = cvdata[:,5]
        self.assertLess(np.max(np.abs(cv - d)), 0.05)
        
if __name__ == "__main__":
    unittest.main()
