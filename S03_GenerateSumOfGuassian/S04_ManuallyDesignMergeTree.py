from S03_Blend2MultivarianceGuassianWithScalarfield import *

if __name__ == '__main__':
    outFolder = join('output', Path(__file__).stem)
    os.makedirs(outFolder, exist_ok=True)

    # Our 2-dimensional distribution will be over variables X and Y
    N = 200
    X = np.linspace(-2, 2, N)
    Y = np.linspace(-2, 2, N)
    X, Y = np.meshgrid(X, Y)

    # Mean vector and covariance matrix
    mu1 = np.array([-1, 0])
    Sigma1 = np.array([[1., 0], [0, 1.]]) / 5

    mu2 = np.array([0.6, 0])
    Sigma2 = np.array([[0.5, 0], [0, 0.5]]) / 5

    mu3 = np.array([1.5, 0])
    Sigma3 = np.array([[0.5, 0], [0, 0.5]]) / 5


    # Pack X and Y into a single 3-dimensional array
    pos = np.empty(X.shape + (2,))
    pos[:, :, 0] = X
    pos[:, :, 1] = Y

    # The distribution on the variables X, Y packed into pos.

    Z1 = multivariate_gaussian(pos, mu1, Sigma1) * 3
    Z2 = multivariate_gaussian(pos, mu2, Sigma2)
    Z3 = multivariate_gaussian(pos, mu3, Sigma3)

    Z = Z1 + Z2 + Z3

    step = 50
    writeOBj(join(outFolder, 'Tree2.obj'), N, X, Y, Z)

    # Tree 1
    mu1 = np.array([-1, 0])
    Sigma1 = np.array([[1., 0], [0, 1.]]) / 5

    mu2 = np.array([1, 0])
    Sigma2 = np.array([[1., 0], [0, 1.]]) / 5

    Z1 = multivariate_gaussian(pos, mu1, Sigma1) * 3
    Z2 = multivariate_gaussian(pos, mu2, Sigma2) * 3

    Z = Z1 + Z2
    writeOBj(join(outFolder, 'Tree1.obj'), N, X, Y, Z)

    obj2vtkFolder(outFolder, )
