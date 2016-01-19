from wham_utils import logSumSlow, logSumFast 
import unittest
import numpy as np




class TestLogSum(unittest.TestCase):
#    def testLogSum(self):
#        vals = np.random.rand(50) + 0.001
#        ls1 = logSumSlow(vals)
#        ls2 = logSumWeave(vals)
#        print "%g - %g = %g" % (ls1, ls2, ls1-ls2)
#        self.assertTrue( abs(ls1 - ls2) < 1e-12, "logSumFast is different from logSumSlow: %g - %g = %g" % (ls1, ls2, ls1-ls2) )
#    def testLogFort(self):
#        vals = np.random.rand(50) + 0.001
#        ls1 = logSumSlow(vals)
#        ls2 = logSumFort(vals)
#        print "%g - %g = %g" % (ls1, ls2, ls1-ls2)
#        self.assertTrue( abs(ls1 - ls2) < 1e-12, "logSumFast is different from logSumSlow: %g - %g = %g" % (ls1, ls2, ls1-ls2) )
    def testLogFast1(self):
        vals = np.random.uniform(1,1000,50) + 0.001
        ls1 = logSumSlow(vals)
        ls2 = logSumFast(vals)
        print "%g - %g = %g" % (ls1, ls2, ls1-ls2)
        self.assertTrue( abs(ls1 - ls2) < 1e-12, "logSumFast is different from logSumSlow: %g - %g = %g" % (ls1, ls2, ls1-ls2) )


if __name__ == "__main__":
    unittest.main()