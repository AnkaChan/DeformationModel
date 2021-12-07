import pyvista as pv
import numpy as np
import itertools
from scipy import sparse

def writeOBj(outObj, N, X, Y, Z):
    '''
    Write a scalar field to a grid mesh
    '''
    file  = open(outObj, 'w')
    for i, j in itertools.product(range(N), range(N)):
        file.write('v %f %f %f\n' %( X[i, j],  Y[i,j], Z[i,j] ))

    for i, j in itertools.product(range(0, N-1), range(1, N)):
        vId = j + i *N
        file.write('f %d %d %d\n' %(vId, vId+1,  vId+N+1, ))
        file.write('f %d %d %d\n' %(vId, vId+N+1,  vId+N, ))

def flatten2DIndex(i,j, gridSize, major='row'):
    '''
    i: row
    j: col
    '''
    if major == 'col':
        return gridSize[0] * j + i
    elif major == 'row':
        return gridSize[1] * i + j

def matchTreeToGrid(treeNodes, flattenGrid):
    numNodes = treeNodes.shape[0]

    matchedVIds = []
    matchDis = []
    for iNode in range(numNodes):
        diff = flattenGrid - treeNodes[iNode, :]
        dis = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)

        matchedVIds.append(np.argmin(dis))
        matchDis.append(dis[matchedVIds[-1]])
        print("Min distance: ", dis[matchedVIds[-1]])

    return matchedVIds


if __name__ == '__main__':
    inputScalarFieldFile = r'X:\Code\CompTopology\Data\S04_ManuallyDesignMergeTree\vtk\ATree1.vtk'
    inputTreeFile = r'X:\Code\CompTopology\Data\S04_ManuallyDesignMergeTree\MergeTree1.vtk'

    gridSize = [200, 200] # [rows, cols]

    scalarField = pv.PolyData(inputScalarFieldFile)
    print(scalarField.points.shape)

    # product('ABCD', repeat=2)
    # AA AB AC AD BA BB BC BD CA CB CC CD DA DB DC DD

    scalarField.points[flatten2DIndex(100, 100, gridSize), 2] = 1.1
    scalarField.points[flatten2DIndex(100, 99, gridSize), 2] = 1
    scalarField.points[flatten2DIndex(100, 101, gridSize), 2] = 1
    scalarField.points[flatten2DIndex(99, 100, gridSize), 2] = 1
    scalarField.points[flatten2DIndex(101, 100, gridSize), 2] = 1

    scalarField.save('test.ply')

    # matching the merge tree node and the grid

    mergeTree = pv.read(inputTreeFile)
    # print(mergeTree.points)

    matchedVId = matchTreeToGrid(mergeTree.points, scalarField.points)
    print(matchedVId)

    # type of 4 nodes
    # saddle points, local maximum, local maximum, root node

    # try to build the Laplacian matrix
    # there will be 199 x 199 terms
    numViables = scalarField.points.shape[0]

    A = np.zeros((numViables, numViables, ), dtype=np.float32)
    # A = np.zeros((numViables, numViables, ), dtype=np.float64)

    for i, j in itertools.product(range(1,199), range(1,199)):
        # y direction second derivatives
        id_ij = flatten2DIndex(i, j, gridSize)
        id_im1j = flatten2DIndex(i-1, j, gridSize)
        id_ip1j = flatten2DIndex(i+1, j, gridSize)

        A[id_ij, id_ij] += 4
        A[id_im1j, id_im1j] += 1
        A[id_ip1j, id_ip1j] += 1

        A[id_ij, id_ip1j] += -2
        A[id_ip1j, id_ij] += -2

        A[id_ij, id_im1j] += -2
        A[id_im1j, id_ij] += -2

        A[id_ip1j, id_im1j] += 1
        A[id_im1j, id_ip1j] += 1

        # x direction second derivatives
        id_ijm1 = flatten2DIndex(i, j-1, gridSize)
        id_ijp1 = flatten2DIndex(i, j+1, gridSize)

        A[id_ij, id_ij] += 4
        A[id_ijm1, id_ijm1] += 1
        A[id_ijp1, id_ijp1] += 1

        A[id_ij, id_ijp1] += -2
        A[id_ijp1, id_ij] += -2

        A[id_ij, id_ijm1] += -2
        A[id_ijm1, id_ij] += -2

        A[id_ijp1, id_ijm1] += 1
        A[id_ijm1, id_ijp1] += 1

    print("number of non-zero elements in A: ", len(np.where(A)[0]), " of ", A.shape[0] * A.shape[1])
    X = scalarField.points[:, 2]
    laplacian = X.transpose() @ A @ X
    print('Laplacian:', laplacian)

    A = sparse.csc_matrix(A)
    sparse.save_npz('LaplacianMat.npz', A)

    np.save('TreeNodes.npy', mergeTree.points)

    # another way to compute Laplacian
    gridData = scalarField.points[:, 2].reshape(gridSize)
    xLaplacian = np.sum((2*gridData[1:-2, 1:-2] - gridData[0:-3, 1:-2] - gridData[2:-1, 1:-2]) ** 2)
    yLaplacian = np.sum((2*gridData[1:-2, 1:-2] - gridData[1:-2, 0:-3,] - gridData[1:-2, 2:-1]) ** 2)

    # Laplacian: 11.521976
    print('Laplacian 2: ', xLaplacian + yLaplacian)