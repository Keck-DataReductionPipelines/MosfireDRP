'''
Convenience functions for Bspline
'''

import numpy as np
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
import CSU

import unittest

def array_to_tuple(A):

        xs = []
        ys = []
        zs = []

        for x in np.arange(A.shape[0]):
                for y in np.arange(A.shape[1]):
                        xs.append(x)
                        ys.append(y)
                        zs.append(A[x,y])

        return list(map(np.array, (xs, ys, zs)))



class TestBsplineFunctions(unittest.TestCase):

        def setUp(self):
                pass

        def test_array_to_tuple(self):
                sa = self.assertTrue

                A = np.array([[1,2], [3,4]])

                (xs,ys,zs) = array_to_tuple(A)
                sa(zs[0] == 1)
                sa(zs[1] == 2)
                sa(zs[2] == 3)
                sa(zs[3] == 4)

                for cnt in range(len(xs)):
                        i = xs[cnt]
                        j = ys[cnt]
                        sa(A[i, j] == zs[cnt])
                        cnt += 1

if __name__ == '__main__':
        unittest.main()

