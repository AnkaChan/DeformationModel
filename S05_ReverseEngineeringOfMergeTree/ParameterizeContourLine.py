from M02_ScalarFieldInterpolation import *
import preprocessing
import matplotlib.pyplot as plt


def findSquareFromEdge(previousEdge, currentEdge, gridSize):
    """
    Compute the two squares adjacent to the current edge and determine the square to use for this step.
    previousEdge: [vertex0, vertex1] (flattened indices), the edge considered from the previous step
    currentEdge: [vertex0, vertex1] (flattened indices), the edge considered in the current step
    Output: three edges in the square excluding the current edge
    """

    i_0, j_0 = to2DIndex(currentEdge[0], gridSize)
    i_1, j_1 = to2DIndex(currentEdge[1], gridSize)

    assert (i_0 - i_1 == 0 and abs(j_0 - j_1) == 1) or (abs(i_0 - i_1) == 1 and j_0 - j_1 == 0)
    # horizontal edge of a square
    if i_0 - i_1 == 0 and abs(j_0 - j_1) == 1:
        square_0 = [[flatten2DIndex(i_0, j_0, gridSize), flatten2DIndex(i_0+1, j_0, gridSize)],
                    [flatten2DIndex(i_0+1, j_0, gridSize), flatten2DIndex(i_1+1, j_1, gridSize)],
                    [flatten2DIndex(i_1+1, j_1, gridSize), flatten2DIndex(i_1, j_1, gridSize)]]
        square_1 = [[flatten2DIndex(i_0, j_0, gridSize), flatten2DIndex(i_0 - 1, j_0, gridSize)],
                    [flatten2DIndex(i_0 - 1, j_0, gridSize), flatten2DIndex(i_1 - 1, j_1, gridSize)],
                    [flatten2DIndex(i_1 - 1, j_1, gridSize), flatten2DIndex(i_1, j_1, gridSize)]]
    # vertical edge of a square
    elif abs(i_0 - i_1) == 1 and j_0 - j_1 == 0:
        square_0 = [[flatten2DIndex(i_0, j_0, gridSize), flatten2DIndex(i_0, j_0-1, gridSize)],
                    [flatten2DIndex(i_0, j_0-1, gridSize), flatten2DIndex(i_1, j_1-1, gridSize)],
                    [flatten2DIndex(i_1, j_1-1, gridSize), flatten2DIndex(i_1, j_1, gridSize)]]
        square_1 = [[flatten2DIndex(i_0, j_0, gridSize), flatten2DIndex(i_0, j_0 + 1, gridSize)],
                    [flatten2DIndex(i_0, j_0 + 1, gridSize), flatten2DIndex(i_1, j_1 + 1, gridSize)],
                    [flatten2DIndex(i_1, j_1 + 1, gridSize), flatten2DIndex(i_1, j_1, gridSize)]]

    if previousEdge in square_0 or [previousEdge[1], previousEdge[0]] in square_0:
        return square_1
    elif previousEdge in square_1 or [previousEdge[1], previousEdge[0]] in square_1:
        return square_0
    else:
        raise ValueError("The previous edge is not in either square that contains the current edge.")


def checkCurrentSquareIntersections(currentSquare, contourEdges):
    intersections = [edge for edge in currentSquare if edge in contourEdges]
    reverseSquare = [[edge[1], edge[0]] for edge in currentSquare]
    reverseIntersections = [edge for edge in reverseSquare if edge in contourEdges]
    if len(intersections) + len(reverseIntersections) > 1:
        raise ValueError("Too many intersection points in the current square.")
    elif len(intersections) + len(reverseIntersections) == 0:
        raise ValueError("No intersection point in the current square.")
    else:
        return True


def isEdge(edge, gridSize):
    i_0, j_0 = to2DIndex(edge[0], gridSize)
    i_1, j_1 = to2DIndex(edge[1], gridSize)

    if check2DCoordValidility(i_0, j_0, gridSize) and check2DCoordValidility(i_1, j_1, gridSize):
        if (i_0 - i_1 == 0 and abs(j_0 - j_1) == 1) or (abs(i_0 - i_1) == 1 and j_0 - j_1 == 0):
            return True
        else:
            return False
    else:
        raise ValueError("At least one of the vertices is not valid.")


def findSaddleNeighborhood(saddleVert, gridSize):
    """
    Find the neighboring edges of a saddle point
    Output: list of 12 edges around a saddle
    """
    i, j = to2DIndex(saddleVert, gridSize)

    left = flatten2DIndex(i, j - 1, gridSize)
    topleft = flatten2DIndex(i + 1, j - 1, gridSize)
    top = flatten2DIndex(i + 1, j, gridSize)
    topright = flatten2DIndex(i + 1, j + 1, gridSize)
    right = flatten2DIndex(i, j + 1, gridSize)
    bottomright = flatten2DIndex(i - 1, j + 1, gridSize)
    bottom = flatten2DIndex(i - 1, j, gridSize)
    bottomleft = flatten2DIndex(i - 1, j - 1, gridSize)

    saddleNeighborEdges = [[left, topleft], [topleft, top], [top, topright], [topright, right], [right, bottomright],
                           [bottomright, bottom], [bottom, bottomleft], [bottomleft, left]]

    return saddleNeighborEdges


def findStartingEdges(saddleVert, saddleNeighborEdges, gridSize, contourEdges, edgesToRemove=None):
    """
    Find the starting currentEdge and starting previousEdge for one saddle
    saddlePoint: a vertex (int) on the grid
    saddleNeighborEdges: 8 edges in the neighborhood of a saddle point but not incident to the saddle point
    edgesToRemove: the edges in the neighborhood that are already starting or ending of the a loop, remove these
    Output: currentEdge - the current edge to start with, previousEdge - the edge used to initialize the process
    """

    if edgesToRemove is not None:
        saddleNeighborEdges = [edge for edge in saddleNeighborEdges if edge not in edgesToRemove]
        reverseEdgesToRemove = [[edge[1], edge[0]] for edge in edgesToRemove]
        saddleNeighborEdges = [edge for edge in saddleNeighborEdges if edge not in reverseEdgesToRemove]

    # determine the first edge to start with
    for iSaddleEdge in saddleNeighborEdges:
        currentEdge = None
        iSaddleEdge_0 = to2DIndex(iSaddleEdge[0], gridSize)
        iSaddleEdge_1 = to2DIndex(iSaddleEdge[1], gridSize)
        if check2DCoordValidility(iSaddleEdge_0[0], iSaddleEdge_0[1], gridSize) and check2DCoordValidility(iSaddleEdge_1[0],iSaddleEdge_1[1], gridSize):
            if iSaddleEdge in contourEdges:
                currentEdge = iSaddleEdge
                break
            elif [iSaddleEdge[1], iSaddleEdge[0]] in contourEdges:
                currentEdge = [iSaddleEdge[1], iSaddleEdge[0]]
                break

    print("initial current edge: ", currentEdge)
    # if there are no more intersections in the saddle's neighborhood, move on to the next saddle
    if currentEdge is None:
        return None, None, None

    currentEdgeIndex = contourEdges.index(currentEdge)

    if isEdge([saddleVert, currentEdge[0]], gridSize):
        previousEdge = [saddleVert, currentEdge[0]]
    elif isEdge([saddleVert, currentEdge[1]], gridSize):
        previousEdge = [saddleVert, currentEdge[1]]
    else:
        raise ValueError("The starting edge is not in the neighborhood of the saddle point.")

    return previousEdge, currentEdge, currentEdgeIndex


def plotEdges(edgesToPlot, step=None):
    """
    Plot the edges up to a specified step.
    """
    x = []
    y = []

    for i, edge in enumerate(edgesToPlot):
        if step is not None and i == step:
            break
        i_0, j_0 = to2DIndex(edge[0], gridSize)
        i_1, j_1 = to2DIndex(edge[1], gridSize)
        x.append(j_0)
        x.append(j_1)
        y.append(i_0)
        y.append(i_1)
    plt.scatter(x, y)
    plt.show()


def reorderContourPointsOneLoop(saddleNeighborEdges, previousEdge, currentEdge, currentEdgeIndex, gridSize, contourEdges, contourWeights, contourHeights):
    """
    Reorder the contour line points, weights and heights for one loop around a saddle point (a saddle can have multiple loops).
    Output: reordered contourEdges, contourWeights, and contourHeights
    """

    if previousEdge == currentEdge or [previousEdge[1], previousEdge[0]] == currentEdge:
        raise ValueError("Previous edge and current edge are the same edge.")

    contourEdgesReordered = [currentEdge]
    contourWeightsReordered = [contourWeights[currentEdgeIndex]]
    contourHeightsReordered = [contourHeights[currentEdgeIndex]]

    # plotEdges([currentEdge, previousEdge])
    currentSquare = findSquareFromEdge(previousEdge, currentEdge, gridSize)

    # print(checkCurrentSquareIntersections(currentSquare, contourEdges))
    if checkCurrentSquareIntersections(currentSquare, contourEdges):
        for edge in currentSquare:
            if edge in contourEdges:
                previousEdge = currentEdge
                currentEdge = edge
                currentEdgeIndex = contourEdges.index(currentEdge)
            elif [edge[1], edge[0]] in contourEdges:
                previousEdge = currentEdge
                currentEdge = [edge[1], edge[0]]
                currentEdgeIndex = contourEdges.index(currentEdge)
    contourEdgesReordered.append(currentEdge)
    contourWeightsReordered.append(contourWeights[currentEdgeIndex])
    contourHeightsReordered.append(contourHeights[currentEdgeIndex])

    while currentEdge not in saddleNeighborEdges and [currentEdge[1], currentEdge[0]] not in saddleNeighborEdges:
        # print("in while loop.")
        currentSquare = findSquareFromEdge(previousEdge, currentEdge, gridSize)
        # print(checkCurrentSquareIntersections(currentSquare, contourEdges))
        if checkCurrentSquareIntersections(currentSquare, contourEdges):
            for edge in currentSquare:
                if edge in contourEdges:
                    previousEdge = currentEdge
                    currentEdge = edge
                    currentEdgeIndex = contourEdges.index(currentEdge)
                elif [edge[1], edge[0]] in contourEdges:
                    previousEdge = currentEdge
                    currentEdge = [edge[1], edge[0]]
                    currentEdgeIndex = contourEdges.index(currentEdge)
        contourEdgesReordered.append(currentEdge)
        contourWeightsReordered.append(contourWeights[currentEdgeIndex])
        contourHeightsReordered.append(contourHeights[currentEdgeIndex])

    return contourEdgesReordered, contourWeightsReordered, contourHeightsReordered


if __name__ == '__main__':
    inputScalarField = './Data/Geodesic.py/larger_domain/MultiGaussian0/MultiGaussian0.obj'
    inputSegmentedField = './Data/Geodesic.py/larger_domain/MultiGaussian0/Seg0.vtk'
    inputMergeTreeNodesField = './Data/Geodesic.py/larger_domain/MultiGaussian0/Node0.vtk'

    gridSize = (256, 257)

    fieldSeg = pv.read(inputSegmentedField)
    fieldSeg.points[:, 2] = fieldSeg['Height']
    field = pv.PolyData(inputScalarField)
    nodes = pv.read(inputMergeTreeNodesField)

    contourInfoPath = './Data/Geodesic.py/larger_domain/MultiGaussian0/ContourInfo'
    contourEdges = preprocessing.readList(os.path.join(contourInfoPath, "contourEdges.txt"))
    contourWeights = preprocessing.readList(os.path.join(contourInfoPath, "contourWeights.txt"))
    contourHeights = preprocessing.readList(os.path.join(contourInfoPath, "contourHeights.txt"))

    # match the merge tree to grid
    iMeshVerts = matchTreeToGrid(nodes.points, fieldSeg.points)

    for iNode in range(nodes.points.shape[0]):
        # skip the nodes that are not saddle points
        print("Saddle point id: ", iNode)
        if nodes['CriticalType'][iNode] != 2:
            continue

        iMeshVert = iMeshVerts[iNode]
        edgesToRemove = []
        saddleNeighborEdges = findSaddleNeighborhood(iMeshVert, gridSize)
        initPEdge, initCEdge, currentEdgeIndex = findStartingEdges(iMeshVert, saddleNeighborEdges, gridSize, contourEdges)

        while initCEdge is not None:
            edgesToRemove.append(initCEdge)
            newEdges, newWeights, newHeights = reorderContourPointsOneLoop(saddleNeighborEdges, initPEdge, initCEdge, currentEdgeIndex, gridSize, contourEdges, contourWeights, contourHeights)
            edgesToRemove.append(newEdges[-1])
            initPEdge, initCEdge, currentEdgeIndex = findStartingEdges(iMeshVert, saddleNeighborEdges, gridSize,
                                                                       contourEdges, edgesToRemove)

            # plot edges up to a certain step number
            # use this to check that the edges are in the correct order and not randomly jumping around
            plotEdges(newEdges, 50)
    # currentEdge = [flatten2DIndex(20, 10, gridSize), flatten2DIndex(20, 11, gridSize)]
    # previousEdge = [flatten2DIndex(18, 11, gridSize), flatten2DIndex(20000, 10000, gridSize)]
    # print(isEdge(currentEdge, gridSize))
    # print(isEdge(previousEdge, gridSize))

    # square = findSquareFromEdge(previousEdge, currentEdge, gridSize)
    # print(square)


