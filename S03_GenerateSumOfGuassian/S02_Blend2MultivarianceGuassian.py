import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import itertools
from os.path import join
from SkelFit.Visualization import obj2vtkFolder

import pyvista as pv

def multivariate_gaussian(pos, mu, Sigma):
    """Return the multivariate Gaussian distribution on array pos."""

    n = mu.shape[0]
    Sigma_det = np.linalg.det(Sigma)
    Sigma_inv = np.linalg.inv(Sigma)
    N = np.sqrt((2*np.pi)**n * Sigma_det)
    # This einsum call calculates (x-mu)T.Sigma-1.(x-mu) in a vectorized
    # way across all the input variables.
    fac = np.einsum('...k,kl,...l->...', pos-mu, Sigma_inv, pos-mu)

    return np.exp(-fac / 2) / N

def writeOBj(outObj, N, X, Y, Z):
    file  = open(outObj, 'w')
    for i, j in itertools.product(range(N), range(N)):
        file.write('v %f %f %f\n' %( X[i, j],  Y[i,j], Z[i,j] ))

    for i, j in itertools.product(range(0, N-1), range(1, N)):
        vId = j + i *N
        file.write('f %d %d %d\n' %(vId, vId+1,  vId+N+1, ))
        file.write('f %d %d %d\n' %(vId, vId+N+1,  vId+N, ))





if __name__ == '__main__':
    outFolder = r'F:\WorkingCopy2\2021_02_15_TestIsometricDeformation\Input\BlendGuassians'
    outFolderVTK = r'F:\WorkingCopy2\2021_02_15_TestIsometricDeformation\Input\BlendGuassiansVTK'
    # Our 2-dimensional distribution will be over variables X and Y
    N = 20
    X = np.linspace(-2, 2, N)
    Y = np.linspace(-2, 2, N)
    X, Y = np.meshgrid(X, Y)

    # Mean vector and covariance matrix
    mu = np.array([0.5, 0.3])
    Sigma = np.array([[1., -0.5], [-0.5, 1.]]) / 3

    mu2 = np.array([-0.3, -0.6])
    Sigma2 = np.array([[0.7, -0.5], [-0.5, 0.7]]) / 3

    # Pack X and Y into a single 3-dimensional array
    pos = np.empty(X.shape + (2,))
    pos[:, :, 0] = X
    pos[:, :, 1] = Y

    # The distribution on the variables X, Y packed into pos.

    zInitial = multivariate_gaussian(pos, np.array([0.0, 0.0]),  np.array([[1., 0], [0, 1.]]) / 3) *3

    Z = multivariate_gaussian(pos, mu, Sigma)

    Z2 = multivariate_gaussian(pos, mu2, Sigma2)

    Z = (Z+Z2) *3

    plt.show()
    # plt.waitkey()

    step = 50

    for i in range(step+1):
        blendedZ = (step-i)/step * zInitial + (i)/step * Z
        writeOBj(join(outFolder, 'A' + str(i).zfill(3) + '.obj',), N, X, Y, blendedZ)

    obj2vtkFolder(outFolder, )

