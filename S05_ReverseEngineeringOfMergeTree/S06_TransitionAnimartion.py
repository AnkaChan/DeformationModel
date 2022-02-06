from S03_SolveScalarFieldContourLIneConstraintInterpolationConstraint import  *
import vtk

def gridCoordToPlaneCoord(i, j, gridSize, xLow=-2, xHigh=2, yLow=-2, yHigh=2):
    return [xLow + (xHigh - xLow) * (i/gridSize[0]), yLow + (yHigh - yLow) * (j/gridSize[1])]


def genScalarField(S1, O2, O3, S2Height, outFile, dType = np.float64):
    O1 = np.array([100, 50])

    O4 = np.array([100, 130])
    S2 = O4

    gridSize = (200, 200)
    criticalPoints = [
        O1,
        O2,
        O3,
        # S1,
        # S2
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
        0,
        0,
        0,
        # 2,
        # 2
    ]
    # contourLineHeight = 2
    # criticalPointHeight = [3, 3, contourLineHeight]

    # Draw the tree
    polyData = vtk.vtkPolyData()
    ptsVtk = vtk.vtkPoints()
    pts = [
        gridCoordToPlaneCoord(O1[1], O1[0], gridSize) + [criticalPointHeight[0]],
        gridCoordToPlaneCoord(O2[1], O2[0], gridSize) + [criticalPointHeight[1]],
        gridCoordToPlaneCoord(O3[1], O3[0], gridSize) + [criticalPointHeight[2]],
        gridCoordToPlaneCoord(S1[1], S1[0], gridSize) + [1],
        gridCoordToPlaneCoord(S2[1], S2[0], gridSize) + [S2Height]
    ]
    for i in range(len(pts)):
        ptsVtk.InsertNextPoint(pts[i])
    polyData.SetPoints(ptsVtk)
    lines = vtk.vtkCellArray()

    edges = [[0, 3], [3, 4], [1, 4], [2, 4]]
    for i in range(len(edges)):
        line = vtk.vtkLine()

        line.GetPointIds().SetId(0, edges[i][0])  # the second 0 is the index of the Origin in the vtkPoints
        line.GetPointIds().SetId(1, edges[i][1])  # the second 1 is the index of P0 in the vtkPoints
        lines.InsertNextCell(line)

    polyData.SetLines(lines)
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polyData)
    writer.SetFileName(join(outFile + '.tree.vtk'))
    writer.Update()

    r1 = np.linalg.norm(O1 - S1)
    r2 = np.linalg.norm(O2 - S2)
    r3 = np.linalg.norm(O3 - S2)
    r4 = np.linalg.norm(O4 - S1)

    coords = []
    for i,j in itertools.product(range(200), range(200)):
        coords .append(np.array([i,j]))
    coords = np.array(coords)

    constraintIds1, constraintWeights1 =  contourLineInterpolationConstraints(O1, r1, gridSize)
    constraintIds4, constraintWeights4 = contourLineInterpolationConstraints(O4, r4, gridSize)

    if  np.linalg.norm(O2 - O3) <= 4:
        constraintIds2 = []
        constraintWeights2 = []
        constraintIds3 = []
        constraintWeights3 = []
    else:
        constraintIds2, constraintWeights2 = contourLineInterpolationConstraints(O2, r2, gridSize)
        constraintIds3, constraintWeights3 = contourLineInterpolationConstraints(O3, r3, gridSize)

    contourLineConstraintIds = constraintIds1 + constraintIds2 + constraintIds3 + constraintIds4
    constraintWeights = constraintWeights1 + constraintWeights2 + constraintWeights3 + constraintWeights4
    constrainHeights = [1 for i in constraintIds1] + [S2Height for i in constraintIds2] + [S2Height for i in constraintIds3] + [1 for i in constraintIds4]

    P = sparse.load_npz('LaplacianMat.npz').astype(dType)
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
        ij = to2DIndex(vId, gridSize)
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
        boundaryVerts.append(flatten2DIndex(i, gridSize[1] - 1, gridSize))

    for j in range(1, gridSize[1] - 1):
        boundaryVerts.append(flatten2DIndex(0, j, gridSize))
        boundaryVerts.append(flatten2DIndex(gridSize[0] - 1, j, gridSize))

    print("Number of boundary constraints", len(boundaryVerts))

    # equality constraints
    equalityConstraints = boundaryVerts + criticalPoints1D
    equalityVal = [0 for i in range(len(boundaryVerts))] + criticalPointHeight

    A = np.zeros((len(equalityConstraints) + len(contourLineConstraintIds), numberOfViables)).astype(dType)
    for i, constraintVId in enumerate(equalityConstraints):
        A[i, constraintVId] = 1

    for i, iConstraint in enumerate(range(len(equalityConstraints), A.shape[0])):
        A[iConstraint, int(contourLineConstraintIds[i][0])] = constraintWeights[i][0]
        A[iConstraint, int(contourLineConstraintIds[i][1])] = constraintWeights[i][1]
        equalityVal.append(constrainHeights[i])

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
    writeOBj(outFile, 200, X, Y, Z)

if __name__ == '__main__':
    # ts = [0.0]
    ts = np.linspace(0,1,50)
    outFolder = join('Data', os.path.basename(__file__))
    os.makedirs(outFolder, exist_ok=True)

    for iStep, t in enumerate(ts):
        S1 = np.array([100, int(70 + t*20)])
        O2 = np.array([100, int(110 + t*20)])
        O3 = np.array([100, int(150 - t*20)])
        S2Height = 1.5 + 0.5*t

        outFile = join(outFolder, 'Animation_' + str(iStep).zfill(3) + '.obj')

        genScalarField(S1, O2, O3, S2Height, outFile)

    obj2vtkFolder(outFolder)