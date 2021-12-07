import pyvista as pv
from S03_SolveScalarFieldContourLIneConstraintInterpolationConstraint import *
from os.path import join
from S02_SolveScalarField import *

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

def check2DCoordValidility(i, j, gridSize):

    return i >= 0 and i < gridSize[0] and j >= 0 and j < gridSize[1]

def findContourLineConstraints(fieldSeg, gridSize, contourLineHeight, edges=None):
    # don't consider boundary edges
    directionsForEdges = [
        (1, 0),
        (0, 1),
    ]
    if edges is None:
        edges = []
        for i,j in tqdm.tqdm(itertools.product(range(gridSize[0]), range(gridSize[1]))):
            for d in directionsForEdges:
                if i + d[0] >= 0 and i + d[0] < gridSize[0] and j + d[1] >= 0 and j + d[1] < gridSize[1]:
                    edges.append([flatten2DIndex(i, j, gridSize=gridSize), flatten2DIndex(i+d[0], j+d[1], gridSize=gridSize)])

    contourEdges = []
    contourLineConstraintWeight = []
    contourLineConstraintHeight = []
    for i, edge in tqdm.tqdm(enumerate(edges)):
        if fieldSeg["SegmentationId"][edge[0]] != fieldSeg["SegmentationId"][edge[1]]:

            contourEdges.append(edge)
            h1 = fieldSeg["Height"][edge[0]]
            h2 = fieldSeg["Height"][edge[1]]
            if h1 > h2:
                if contourLineHeight.get(
                        (fieldSeg["SegmentationId"][edge[0]], fieldSeg["SegmentationId"][edge[1]])) is None:
                    continue
                contourHeight = contourLineHeight[fieldSeg["SegmentationId"][edge[0]], fieldSeg["SegmentationId"][edge[1]]]
            else:
                if contourLineHeight.get(
                        (fieldSeg["SegmentationId"][edge[1]], fieldSeg["SegmentationId"][edge[0]])) is None:
                    continue
                contourHeight = contourLineHeight[
                    fieldSeg["SegmentationId"][edge[1]], fieldSeg["SegmentationId"][edge[0]]]

            w1 = (contourHeight - h2) / (h1 - h2)
            assert 0 <= w1 <=1
            w2 = 1 - w1
            contourLineConstraintWeight.append([w1, w2])
            contourLineConstraintHeight.append(contourHeight)

    return contourEdges, contourLineConstraintWeight, contourLineConstraintHeight

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


def findContourLineHeight(nodes, fieldSeg, gridSize):
    contourLineHeight = {} # key (a, b): a: higher segment, c lower segment
    iMeshVerts = matchTreeToGrid(nodes.points, fieldSeg.points)

    for iNode in range(nodes.points.shape[0]):

        if nodes['CriticalType'][iNode] != 2:
            continue
        # how to find the height of a contour line?
        # the height of a contour line is the height of a saddle point where to contour line meets
        # at the neighborhood of each saddle point, there will be three types of nodes: two higher segment and one lower segment
        # we will need to find the two higher seg, say (a,b) and the lower seg c
        # the height of contour line a-c, b-c will be the height (scalar) of this saddle point

        # match node to grid mesh
        iMeshVert = iMeshVerts[iNode]
        i, j = to2DIndex(iMeshVert, gridSize)

        segsInNeighborhood = []
        heights = []

        for d in directions:
            if check2DCoordValidility(i+d[0], j+d[1], gridSize):
                neiVertId = flatten2DIndex(i+d[0], j+d[1], gridSize)
                if fieldSeg['SegmentationId'][neiVertId] not in segsInNeighborhood:
                    segsInNeighborhood.append(fieldSeg['SegmentationId'][neiVertId])
                    heights.append(fieldSeg['Height'][neiVertId])

        assert len(segsInNeighborhood) == 3
        sortedId = np.argsort(heights)
        contourLineHeight[(segsInNeighborhood[sortedId[2]], segsInNeighborhood[sortedId[0]])] = nodes['Scalar'][iNode]
        contourLineHeight[(segsInNeighborhood[sortedId[1]], segsInNeighborhood[sortedId[0]])] = nodes['Scalar'][iNode]

    return contourLineHeight


if __name__ == '__main__':
    inputScalarField = r'X:\Code\CompTopology\Data\S04_GenerateNewScalarField.py\M00000.obj'
    inputSegmentedField = r'X:\Code\CompTopology\Data\S04_GenerateNewScalarField.py\Topology\Seg.vtk'
    inputMergeTreeNodesField = r'X:\Code\CompTopology\Data\S04_GenerateNewScalarField.py\Topology\Node.vtk'
    gridSize = (200, 200)
    dType = np.float64
    P = sparse.load_npz('LaplacianMat.npz').astype(dType)

    outFolder = r'./Data/' + os.path.basename(__file__)
    os.makedirs(outFolder, exist_ok=True)
    numberOfVariable = gridSize[0] * gridSize[1]

    fieldSeg = pv.read(inputSegmentedField)
    fieldSeg.points[:, 2] = fieldSeg['Height']
    fieldSeg.save(join(outFolder, 'ScalarFieldSimplified.vtk'))
    field = pv.PolyData(inputScalarField)

    nodes = pv.read(inputMergeTreeNodesField)
    contourLineHeight = findContourLineHeight(nodes, fieldSeg, gridSize)

    # print(fieldSeg.array_names)
    # print(fieldSeg.field_arrays["SegmentationId"])
    # fieldSeg["SegmentationId"]

    # find the contour line constraint
    # by iterating all the edges to find out the edge that cross two segmentations

    contourEdges, contourLineConstraintWeight, contourLineConstraintHeight = findContourLineConstraints(fieldSeg, gridSize, contourLineHeight)

    contourLinePoints = []

    for edge, weights in zip(contourEdges, contourLineConstraintWeight):
        contourLinePoints.append(fieldSeg.points[edge[0],:]*weights[0] + fieldSeg.points[edge[1],:]*weights[1])

    contourLinePointsCloud = pv.PolyData()
    contourLinePointsCloud.points = np.array(contourLinePoints)

    contourLinePointsCloud.save('contourLine.ply')

    # inequality constraints
    ## match critical points to grid verts
    iMeshVerts = matchTreeToGrid(nodes.points, fieldSeg.points)
    ## critical points
    G = []
    h = []
    inequilityConstraintIds = []
    for i, vId in enumerate(iMeshVerts):
        ij = to2DIndex(vId,gridSize)
        # print(flatten2DIndex(ij[0], ij[1],gridSize))
        criticalType = nodes['CriticalType'][i]
        if criticalType == 3:
            for d in directions:
                row = np.zeros((numberOfVariable,))

                neighborIj = (ij[0] + d[0], ij[1] + d[1])
                neighborVId = flatten2DIndex(neighborIj[0], neighborIj[1], gridSize)
                inequilityConstraintIds.append(neighborVId)

                row[vId] = -1
                row[neighborVId] = 1

                G.append(row)
                h.append(0)
    G = np.array(G).astype(dType)
    h = np.array(h).astype(dType)

    # equality constraints
    ## boundary vertices
    boundaryVerts = []
    boundaryVertsHeight = []
    for i in range(gridSize[0]):
        boundaryVerts.append(flatten2DIndex(i, 0, gridSize))
        boundaryVertsHeight.append(fieldSeg["Height"][flatten2DIndex(i, 0, gridSize)])

        boundaryVerts.append(flatten2DIndex(i, gridSize[1]-1, gridSize))
        boundaryVertsHeight.append(fieldSeg["Height"][flatten2DIndex(i, gridSize[1]-1, gridSize)])

    for j in range(1, gridSize[1]-1):
        boundaryVerts.append(flatten2DIndex(0, j, gridSize))
        boundaryVertsHeight.append(fieldSeg["Height"][flatten2DIndex(0, j, gridSize)])

        boundaryVerts.append(flatten2DIndex(gridSize[0]-1, j, gridSize))
        boundaryVertsHeight.append(fieldSeg["Height"][flatten2DIndex(gridSize[0]-1, j, gridSize)])


    criticalPoints1D = []
    criticalPointHeight = []
    for i, iVert in enumerate(iMeshVerts):
        if nodes['CriticalType'][i] == 3 or nodes['CriticalType'][i] == 2:
            criticalPoints1D.append(iVert)
            criticalPointHeight.append(nodes['Scalar'][i])

    equalityConstraints = boundaryVerts + criticalPoints1D
    equalityVal = boundaryVertsHeight + criticalPointHeight

    A = np.zeros((len(equalityConstraints) + len(contourEdges), numberOfVariable)).astype(dType)
    for i, constraintVId in enumerate(equalityConstraints):
        A[i, constraintVId] = 1

    for i, iConstraint in enumerate(range(len(equalityConstraints), A.shape[0])):

        A[iConstraint, int(contourEdges[i][0])] = contourLineConstraintWeight[i][0]
        A[iConstraint, int(contourEdges[i][1])] = contourLineConstraintWeight[i][1]
        equalityVal.append(contourLineConstraintHeight[i])

    b = np.array(equalityVal).astype(dType)
    q = np.zeros((numberOfVariable,)).astype(dType)

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
    X = np.linspace(-1, 1, N)
    Y = np.linspace(-1, 1, N)
    X, Y = np.meshgrid(X, Y)
    writeOBj(join(outFolder, 'result.obj'), 200, X, Y, Z)
    obj2vtkFolder(outFolder)
