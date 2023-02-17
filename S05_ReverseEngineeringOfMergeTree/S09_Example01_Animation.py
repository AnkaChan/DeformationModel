'''
This example interpolate a transition animation guided by a merge tree animation.
'''

from M01_TopologicalExtraction import *
from M02_ScalarFieldInterpolation import *

def gridCoordToPlaneCoord(i, j, gridSize, xLow=-2, xHigh=2, yLow=-2, yHigh=2):
    '''
        convert grid coodinates (i, j) to actual (x, y) coordinates on the plane
    '''
    return [xLow + (xHigh - xLow) * (i/gridSize[0]), yLow + (yHigh - yLow) * (j/gridSize[1])]

def makeMergeTree(S1, O2, O3, S2Height, gridSize):
    O1 = np.array([100, 50]) # this leaf node is fixed
    O4 = np.array([100, 130]) # O2 and O3 will merge to this node if they are close enough
    S2 = O4
    edges = [[0, 3], [3, 4], [1, 4], [2, 4]] # edges for this tree

    criticalPoints = [
        O1,
        O2,
        O3,
    ]
    criticalPointHeight = [
        2,
        2,
        2,
        # 1,
        # 1.5
    ]
    criticalPoints1D = [flatten2DIndex(coord[0], coord[1], gridSize) for coord in criticalPoints]
    # 0 local maximum; 1 local minimum; 2 saddle
    criticalPointType = [
        0, # 0: local maximum
        0,
        0,
        # 2,
        # 2
    ]

    pts = [
        gridCoordToPlaneCoord(O1[1], O1[0], gridSize) + [criticalPointHeight[0]],
        gridCoordToPlaneCoord(O2[1], O2[0], gridSize) + [criticalPointHeight[1]],
        gridCoordToPlaneCoord(O3[1], O3[0], gridSize) + [criticalPointHeight[2]],
        gridCoordToPlaneCoord(S1[1], S1[0], gridSize) + [1],
        gridCoordToPlaneCoord(S2[1], S2[0], gridSize) + [S2Height]
    ]

    drawMergeTree(pts, edges, join(outFile + '.tree.vtk'))

    # the contour line are just circles
    # and their radius are computed automatically based on the position of saddle and leaf points
    r1 = np.linalg.norm(O1 - S1)
    r2 = np.linalg.norm(O2 - S2)
    r3 = np.linalg.norm(O3 - S2)
    r4 = np.linalg.norm(O4 - S1)

    # contour line constraints
    coords = []
    for i,j in itertools.product(range(200), range(200)):
        coords .append(np.array([i,j]))
    coords = np.array(coords)

    constraintIds1, constraintWeights1 =  contourLineInterpolationConstraints(O1, r1, gridSize)
    constraintIds4, constraintWeights4 = contourLineInterpolationConstraints(O4, r4, gridSize)

    if np.linalg.norm(O2 - O3) <= 4:
        # if O2 and O3 are close enough, merge them, and remove the contour line constraints
        constraintIds2 = []
        constraintWeights2 = []
        constraintIds3 = []
        constraintWeights3 = []
    else:
        constraintIds2, constraintWeights2 = contourLineInterpolationConstraints(O2, r2, gridSize)
        constraintIds3, constraintWeights3 = contourLineInterpolationConstraints(O3, r3, gridSize)

    equalityConstraints = constraintIds1 + constraintIds2 + constraintIds3 + constraintIds4
    equalityConstraintWeights = constraintWeights1 + constraintWeights2 + constraintWeights3 + constraintWeights4
    equalityConstraintVals = [1 for i in constraintIds1] + [S2Height for i in constraintIds2] + [S2Height for i in constraintIds3] + [1 for i in constraintIds4]

    # extremity constraints: includes both equality and inequality
    # convert extremity constraints to  equality constraints and inequality constraints

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
    inequalityConstraints = []
    inequalityConstraintWeights = []
    inequalityConstraintVals = []

    for i, vId in enumerate(criticalPoints1D):
        ij = to2DIndex(vId, gridSize)
        # print(flatten2DIndex(ij[0], ij[1],gridSize))
        # equality, value of leaf node
        equalityConstraints.append([vId])
        equalityConstraintWeights.append([1.])
        equalityConstraintVals.append(criticalPointHeight[i])

        criticalType = criticalPointType[i]
        if criticalType == 0:
            for d in directions:
                neighborIj = (ij[0] + d[0], ij[1] + d[1])
                neighborVId = flatten2DIndex(neighborIj[0], neighborIj[1], gridSize)

                # inequality, larger than neighbors
                inequalityConstraints.append([vId, neighborVId, ])
                inequalityConstraintWeights.append([-1., 1.])
                inequalityConstraintVals.append(0.)

    return equalityConstraints, equalityConstraintWeights, equalityConstraintVals, inequalityConstraints, inequalityConstraintWeights, inequalityConstraintVals

def genScalarField(P, S1, O2, O3, S2Height, outFile, gridSize, dType = np.float64):
    '''
        O*: position of leaf nodes
        S*: position of saddle points
    '''

    equalityConstraints, equalityConstraintWeights, equalityConstraintVals, inequalityConstraints, inequalityConstraintWeights, inequalityConstraintVals = makeMergeTree(S1, O2, O3, S2Height, gridSize)
    print(len(equalityConstraints))
    print(len(equalityConstraints[500]))
    print(len(equalityConstraints[500]))
    print(len(equalityConstraintWeights))
    print(len(equalityConstraintVals))
    print(equalityConstraints[500])
    print(equalityConstraintWeights[500])
    print(equalityConstraintVals[500])
    Z = interpolateScalarField(gridSize, P, equalityConstraints, equalityConstraintWeights, equalityConstraintVals, inequalityConstraints, inequalityConstraintWeights, inequalityConstraintVals, fixBoundary=True)

    N = gridSize[1]
    X = np.linspace(-2, 2, N)
    Y = np.linspace(-2, 2, N)
    X, Y = np.meshgrid(X, Y)
    # writeOBj(outFile, 200, X, Y, Z)
    writeOBj(outFile, X, Y, Z, gridSize)


if __name__ == '__main__':
    # timeSteps = np.linspace(0,1,50)
    # outFolder = join('Data', os.path.basename(__file__)) + '/'
    # os.makedirs(outFolder, exist_ok=True)
    #
    # gridSize = (200, 200)
    # P = genLaplacianMat(gridSize)
    #
    # for iStep, t in enumerate(timeSteps):
    #     S1 = np.array([100, int(70 + t*20)])
    #     O2 = np.array([100, int(110 + t*20)])
    #     O3 = np.array([100, int(150 - t*20)])
    #     S2Height = 1.5 + 0.5*t
    #
    #     outFile = join(outFolder, 'Animation_' + str(iStep).zfill(3) + '.obj')
    #
    #     genScalarField(P, S1, O2, O3, S2Height, outFile, gridSize=gridSize)
    #     break
    # print(outFolder)
    # obj2vtkFolder(outFolder, outVtkFolder=outFolder)

    print(gridCoordToPlaneCoord(1, 121, (200, 200)))
