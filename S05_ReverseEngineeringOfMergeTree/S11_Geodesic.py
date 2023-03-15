from M02_ScalarFieldInterpolation import *
import M03_Preprocessing

if __name__ == '__main__':
    # obj2vtkFolder('./Data/S11_Geodesic.py/larger_domain/', outVtkFolder='./Data/S11_Geodesic.py/larger_domain/')

    # inputScalarField = './Data/S11_Geodesic.py/larger_domain/MultiGaussian0.obj'
    #
    # inputSegmentedField = './Data/S11_Geodesic.py/larger_domain/Seg0.vtk'
    # inputMergeTreeNodesField = './Data/S11_Geodesic.py/larger_domain/Node0.vtk'
    inputScalarField = './Data/S11_Geodesic.py/larger_domain/MultiGaussian1/MultiGaussian1.obj'
    inputSegmentedField = './Data/S11_Geodesic.py/larger_domain/MultiGaussian1/Seg1.vtk'
    inputMergeTreeNodesField = './Data/S11_Geodesic.py/larger_domain/MultiGaussian1/Node1.vtk'
    inputMergeTreeEdges = './Data/S11_Geodesic.py/larger_domain/MultiGaussian1/Node1.vtk'
    # inputScalarField = './Data/S11_Geodesic.py/larger_domain/MultiGaussian0/MultiGaussian0.obj'
    # inputSegmentedField = './Data/S11_Geodesic.py/larger_domain/MultiGaussian0/Seg0.vtk'
    # inputMergeTreeNodesField = './Data/S11_Geodesic.py/larger_domain/MultiGaussian0/Node0.vtk'
    # inputScalarField = './Data/S11_Geodesic.py/original/MultiGaussian0/MultiGaussian0.obj'
    # inputSegmentedField = './Data/S11_Geodesic.py/original/MultiGaussian0/Seg0.vtk'
    # inputMergeTreeNodesField = './Data/S11_Geodesic.py/original/MultiGaussian0/Node0.vtk'
    # inputScalarField = './Data/S11_Geodesic.py/original/MultiGaussian1/MultiGaussian1.obj'
    # inputSegmentedField = './Data/S11_Geodesic.py/original/MultiGaussian1/Seg1.vtk'
    # inputMergeTreeNodesField = './Data/S11_Geodesic.py/original/MultiGaussian1/Node1.vtk'



    gridSize = (256, 257) #(Y, X) resolution in paraview (rows + 1, cols + 1)
    dType = np.float64
    P = sparse.csc_matrix(genLaplacianMat(gridSize))
    print(P.shape)
    directions = generate_directions(3)

    outFolder = r'./Data/' + os.path.basename(__file__) + '/larger_domain/MultiGaussian1/results/'
    # print(outFolder)
    os.makedirs(outFolder, exist_ok=True)
    numberOfVariable = gridSize[0] * gridSize[1]

    fieldSeg = pv.read(inputSegmentedField)
    fieldSeg.points[:, 2] = fieldSeg['Height']
    # fieldSeg.save(join(outFolder, 'ScalarFieldSimplified.vtk'))
    field = pv.PolyData(inputScalarField)

    nodes = pv.read(inputMergeTreeNodesField)
    contourLineHeight = findContourLineHeight(nodes, fieldSeg, gridSize, directions)
    # print(contourLineHeight)

    contourEdges, contourLineConstraintWeight, contourLineConstraintHeight = findContourLineConstraints(fieldSeg, gridSize,
                                                                                                        contourLineHeight)
    contourInfoPath = './Data/S11_Geodesic.py/larger_domain/MultiGaussian1/ContourInfo'
    os.makedirs(contourInfoPath, exist_ok=True)
    preprocessing.saveList(os.path.join(contourInfoPath, "contourEdges.txt"), contourEdges)
    preprocessing.saveList(os.path.join(contourInfoPath, "contourWeights.txt"), contourLineConstraintWeight)
    preprocessing.saveList(os.path.join(contourInfoPath, "contourHeights.txt"), contourLineConstraintHeight)


    # print(contourEdges)
    # print(contourLineConstraintHeight)
    contourLinePoints = []

    # save and visualize contour line
    for edge, weights in zip(contourEdges, contourLineConstraintWeight):
        contourLinePoints.append(fieldSeg.points[edge[0], :] * weights[0] + fieldSeg.points[edge[1], :] * weights[1])

    contourLinePointsCloud = pv.PolyData()
    contourLinePointsCloud.points = np.array(contourLinePoints)
    contourLinePointsCloud.save(join(outFolder, 'contourLine.ply'))


    # the boundary value constraints, keep boundary value as it is
    boundaryVerts = []
    boundaryVertsWeights = []
    boundaryVertsHeight = []

    for i in range(gridSize[0]):
        boundaryVerts.append([flatten2DIndex(i, 0, gridSize)])
        boundaryVertsWeights.append([1.])
        boundaryVertsHeight.append(fieldSeg["Height"][flatten2DIndex(i, 0, gridSize)])

        boundaryVerts.append([flatten2DIndex(i, gridSize[1] - 1, gridSize)])
        boundaryVertsWeights.append([1.])
        boundaryVertsHeight.append(fieldSeg["Height"][flatten2DIndex(i, gridSize[1] - 1, gridSize)])

    equalityConstraints = contourEdges + boundaryVerts
    equalityConstraintWeights = contourLineConstraintWeight + boundaryVertsWeights
    equalityConstraintVals = contourLineConstraintHeight + boundaryVertsHeight

    for j in range(1, gridSize[1] - 1):
        boundaryVerts.append([flatten2DIndex(0, j, gridSize)])
        boundaryVertsWeights.append([1.])
        boundaryVertsHeight.append(fieldSeg["Height"][flatten2DIndex(0, j, gridSize)])

        boundaryVerts.append([flatten2DIndex(gridSize[0] - 1, j, gridSize)])
        boundaryVertsWeights.append([1.])
        boundaryVertsHeight.append(fieldSeg["Height"][flatten2DIndex(gridSize[0] - 1, j, gridSize)])

    # inequality constraints
    ## match critical points to grid verts
    iMeshVerts = matchTreeToGrid(nodes.points, fieldSeg.points)
    ## critical points
    inequalityConstraints = []
    inequalityConstraintWeights = []
    inequalityConstraintVals = []
    for i, vId in enumerate(iMeshVerts):
        ij = to2DIndex(vId, gridSize)
        # print(flatten2DIndex(ij[0], ij[1],gridSize))
        criticalType = nodes['CriticalType'][i]
        # equality constraints for extremity
        equalityConstraints.append([vId])
        equalityConstraintWeights.append([1.])
        equalityConstraintVals.append(nodes['Scalar'][i])
        if criticalType == 3:
            # 3: merge tree, local maximum
            for d in directions:
                neighborIj = (ij[0] + d[0], ij[1] + d[1])
                neighborVId = flatten2DIndex(neighborIj[0], neighborIj[1], gridSize)

                # inequality, larger than neighbors
                inequalityConstraints.append([vId, neighborVId, ])
                inequalityConstraintWeights.append([-1., 1.])
                inequalityConstraintVals.append(0.)

    Z = interpolateScalarField(gridSize, P, equalityConstraints, equalityConstraintWeights, equalityConstraintVals,
                               inequalityConstraints, inequalityConstraintWeights, inequalityConstraintVals,
                               fixBoundary=False)
    Nx = gridSize[0]
    Ny = gridSize[1]

    X = np.linspace(min(field.points[:, 0]), max(field.points[:, 0]), Ny)
    Y = np.linspace(min(field.points[:, 1]), max(field.points[:, 1]), Nx)

    # np.meshgrid default indexing='xy' treat X[j,i], Y[j,i]
    X, Y = np.meshgrid(X, Y)
    writeOBj(join(outFolder, 'result.obj'), X, Y, Z, gridSize)