# from S01_PrepareData import *
import numpy as np
import itertools, glob, os, tqdm
from pathlib import Path
import pyvista as pv
from os.path import join

from qpsolvers import solve_qp, available_solvers
from scipy import sparse
print('Avaliable qp solvers: ', available_solvers)
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

def to2DIndex(id, gridSize, major='row'):
    '''
    i: row
    j: col
    '''
    if major == 'col':
        return (int(id % gridSize[0]), int((np.floor(id / gridSize[0]))))
    elif major == 'row':
        return  (int(np.floor(id / gridSize[1])), int(id % gridSize[1]))
def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

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

def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)

def contourLineConstraint(O, r, disThreshold=0.5):
    coords = []
    for i,j in itertools.product(range(200), range(200)):
        if i <= O[0] + r and i >= O[0] - r and j <= O[1] + r and j >= O[1] - r:
            coords .append(np.array([i,j]))
    coords = np.array(coords)
    # determine if pixel's distance is less than 0.5
    ## x direction
    toCircleDisXDir = np.abs(np.sqrt(r**2 - (O[1] - coords[:, 1])**2) - np.abs(O[0] - coords[:, 0]))
    toCircleDisYDir = np.abs(np.sqrt(r**2 - (O[0] - coords[:, 0])**2) - np.abs(O[1] - coords[:, 1]))

    pixelsOnCircle = np.union1d(np.where(toCircleDisXDir<disThreshold)[0], np.where(toCircleDisYDir<disThreshold)[0])
    coordsOnCirlcle = coords[pixelsOnCircle]

    return coordsOnCirlcle

def contourLineConstraintOnInputCoords(coords, O, r, disThreshold=0.5):
    # determine if pixel's distance is less than 0.5
    ## x direction
    toCircleDisXDir = np.abs(np.sqrt(r**2 - (O[1] - coords[:, 1])**2) - np.abs(O[0] - coords[:, 0]))
    toCircleDisYDir = np.abs(np.sqrt(r**2 - (O[0] - coords[:, 0])**2) - np.abs(O[1] - coords[:, 1]))

    pixelsOnCircle = np.union1d(np.where(toCircleDisXDir<disThreshold)[0], np.where(toCircleDisYDir<disThreshold)[0])

    return pixelsOnCircle

def contourLineInterpolationConstraints(O, r, gridSize = (200, 200)):

    constraintIds = []
    constraintWeights = []

    # intersection with axis 0
    for i in range(int(np.ceil(O[0]-r)), int(np.floor(O[0]+r))):
        dy = np.sqrt(r**2 - (i-O[0])**2)
        if dy > 1e-5:
            j = O[1] + dy
            constraintIds.append([flatten2DIndex(i, np.floor(j), gridSize=gridSize), flatten2DIndex(i, np.ceil(j), gridSize=gridSize)])
            constraintWeights.append([1-(j-np.floor(j)), 1-(np.ceil(j)-j)])

            j = O[1] - dy
            constraintIds.append([flatten2DIndex(i, np.floor(j), gridSize=gridSize), flatten2DIndex(i, np.ceil(j), gridSize=gridSize)])
            constraintWeights.append([1-(j-np.floor(j)), 1-(np.ceil(j)-j)])
        else:
            j = O[1]
            constraintIds.append(
                [flatten2DIndex(i, np.floor(j), gridSize=gridSize), flatten2DIndex(i, np.ceil(j), gridSize=gridSize)])
            constraintWeights.append([1 - (j - np.floor(j)), 1 - (np.ceil(j)-j)])

    # intersection with axis 0
    for j in  range(int(np.ceil(O[1]-r)), int(np.floor(O[1]+r))):
        dy = np.sqrt(r ** 2 - (j - O[1]) ** 2)
        if dy > 1e-5:
            i = O[0] + dy
            constraintIds.append([flatten2DIndex(np.floor(i), j, gridSize=gridSize),
                                  flatten2DIndex(np.ceil(i), j, gridSize=gridSize)])
            constraintWeights.append([1 - (i - np.floor(i)), 1 - (np.ceil(i)-i)])

            i = O[0] - dy
            constraintIds.append([flatten2DIndex(np.floor(i), j, gridSize=gridSize),
                                  flatten2DIndex(np.ceil(i), j, gridSize=gridSize)])
            constraintWeights.append([1 - (i - np.floor(i)), 1 - (np.ceil(i)-i)])
        else:
            i = O[0]
            constraintIds.append([flatten2DIndex(np.floor(i), j, gridSize=gridSize),
                                  flatten2DIndex(np.ceil(i), j, gridSize=gridSize)])
            constraintWeights.append([1 - (i - np.floor(i)), 1 - (np.ceil(i) - i)])

    return constraintIds, constraintWeights

def obj2vtkFolder(inObjFolder, inFileExt='obj', outVtkFolder=None, processInterval=[], addFaces=False,
                  addABeforeName=True, faceMesh=None, ):
    # addFaces = True
    if outVtkFolder is None:
        outVtkFolder = inObjFolder + r'\vtk'

    objFiles = glob.glob(inObjFolder + r'\*.' + inFileExt)


    meshWithFaces = None

    os.makedirs(outVtkFolder, exist_ok=True)
    if len(processInterval) == 2:
        objFiles = objFiles[processInterval[0]: processInterval[1]]
    for f in tqdm.tqdm(objFiles, desc=inFileExt + " to vtk"):
        fp = Path(f)

        mesh = pv.read(f)
        if addFaces and meshWithFaces is not None:
            mesh.faces = meshWithFaces.faces
        # else:
        #     mesh.faces = np.empty((0,), dtype=np.int32)
        if addABeforeName:
            outName = outVtkFolder + r'\\A' + fp.stem + '.vtk'
        else:
            outName = outVtkFolder + r'\\' + fp.stem + '.vtk'

        mesh.point_arrays['Height'] = mesh.points[:, 2]

        mesh.save(outName)

if __name__ == '__main__':
    # compute the constraint on the contour line
    O1 = np.array([100, 70])
    O2 = np.array([100, 130])
    saddlePoint = np.array([100,100])
    gridSize = (200, 200)
    criticalPoints = [
        [100, 70],
        [100, 130],
        saddlePoint
    ]
    criticalPoints1D = [flatten2DIndex(coord[0], coord[1], gridSize) for coord in criticalPoints]
    # 0 local maximum; 1 local minimum; 2 saddle
    criticalPointType = [
        0,
        0,
        2
    ]
    contourLineHeight=2
    criticalPointHeight = [3,3,contourLineHeight]
    dType = np.float64

    #
    r1 = np.linalg.norm(O1-saddlePoint)
    r2 = np.linalg.norm(O2-saddlePoint)

    constraintIds1, constraintWeights1 =  contourLineInterpolationConstraints(O1, r1, gridSize)
    constraintIds2, constraintWeights2 =  contourLineInterpolationConstraints(O2, r2, gridSize)

    constraintIds = constraintIds1 + constraintIds2
    constraintWeights = constraintWeights1 + constraintWeights2

    P = sparse.load_npz('LaplacianMat.npz').astype(dType)
    # print('P is symmetric:', check_symmetric(P))

    gridSize = [200, 200]
    numberOfViables = P.shape[0]

    # inequality constraint
    G = []
    h = []
    inequilityConstraintIds = []
    directions = [
        (-1, -1),
        (0, -1),
        (1, -1),
        (-1, 0),
        (1, 0),
        (-1, 1),
        (0, 1),
        (1, 1),
    ]
    for i, vId in enumerate(criticalPoints1D):
        ij = to2DIndex(vId,gridSize)
        # print(flatten2DIndex(ij[0], ij[1],gridSize))
        criticalType = criticalPointType[i]
        if criticalType == 0:
            for d in directions:
                row = np.zeros((numberOfViables,))

                neighborIj = (ij[0] + d[0], ij[1] + d[1])
                neighborVId = flatten2DIndex(neighborIj[0], neighborIj[1], gridSize)
                inequilityConstraintIds.append(neighborVId)

                row[vId] = -1
                row[neighborVId] = 1

                G.append(row)
                h.append(0)

    G = np.array(G).astype(dType)
    h = np.array(h).astype(dType)

    # get boundary vertices
    boundaryVerts = []
    for i in range(gridSize[0]):
        boundaryVerts.append(flatten2DIndex(i, 0, gridSize))
        boundaryVerts.append(flatten2DIndex(i, gridSize[1]-1, gridSize))

    for j in range(1, gridSize[1]-1):
        boundaryVerts.append(flatten2DIndex(0, j, gridSize))
        boundaryVerts.append(flatten2DIndex(gridSize[0]-1, j, gridSize))

    print("Number of boundary constraints", len(boundaryVerts))

    # equality constraints
    equalityConstraints = boundaryVerts + criticalPoints1D
    equalityVal = [0 for i in range(len(boundaryVerts))] + criticalPointHeight

    A = np.zeros((len(equalityConstraints) + len(constraintIds), numberOfViables)).astype(dType)
    for i, constraintVId in enumerate(equalityConstraints):
        A[i, constraintVId] = 1

    for i, iConstraint in enumerate(range(len(equalityConstraints), A.shape[0])):

        A[iConstraint, int(constraintIds[i][0])] = constraintWeights[i][0]
        A[iConstraint, int(constraintIds[i][1])] = constraintWeights[i][1]
        equalityVal.append(contourLineHeight)

    b = np.array(equalityVal).astype(dType)

    q = np.zeros((numberOfViables,)).astype(dType)

    # P = sparse.csc_matrix(P)
    G = sparse.csc_matrix(G)
    # q = sparse.csc_matrix(q)
    # h = sparse.csc_matrix(h)
    A = sparse.csc_matrix(A)
    # b = sparse.csc_matrix(b)

    Z = solve_qp(P, q, G, h, A, b, verbose=True, solver='osqp')
    # x = solve_qp(P, q, G, h, A, b, verbose=True, solver='ecos')

    Z=Z.reshape(gridSize)
    N = 200
    X = np.linspace(-2, 2, N)
    Y = np.linspace(-2, 2, N)
    X, Y = np.meshgrid(X, Y)
    writeOBj('result.obj', 200, X, Y, Z)

    # obj2vtkFolder('.')
    # print(x)