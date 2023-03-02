# from M02_ScalarFieldInterpolation import *
from M01_TopologicalExtraction import *

import preprocessing
import matplotlib.pyplot as plt


if __name__ == '__main__':
    inputScalarField = './Data/S11_Geodesic.py/larger_domain/MultiGaussian0/MultiGaussian0.obj'
    inputSegmentedField = './Data/S11_Geodesic.py/larger_domain/MultiGaussian0/Seg0.vtk'
    inputMergeTreeNodesField = './Data/S11_Geodesic.py/larger_domain/MultiGaussian0/Node0.vtk'

    gridSize = (256, 257)

    fieldSeg = pv.read(inputSegmentedField)
    fieldSeg.points[:, 2] = fieldSeg['Height']
    field = pv.PolyData(inputScalarField)
    nodes = pv.read(inputMergeTreeNodesField)

    contourInfoPath = './Data/S11_Geodesic.py/larger_domain/MultiGaussian0/ContourInfo'
    contourEdges = preprocessing.readList(os.path.join(contourInfoPath, "contourEdges.txt"))
    contourWeights = preprocessing.readList(os.path.join(contourInfoPath, "contourWeights.txt"))
    contourHeights = preprocessing.readList(os.path.join(contourInfoPath, "contourHeights.txt"))

    contourEdges = [(e[0], e[1]) for e in contourEdges]

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

        saddleAllContourEdges = []
        saddleAllContourWeights = []
        saddleAllContourHeights = []

        while initCEdge is not None:
            edgesToRemove.append(initCEdge)
            newEdges, newWeights, newHeights = reorderContourPointsOneLoop(saddleNeighborEdges, initPEdge, initCEdge, currentEdgeIndex, gridSize, contourEdges, contourWeights, contourHeights)
            edgesToRemove.append(newEdges[-1])
            initPEdge, initCEdge, currentEdgeIndex = findStartingEdges(iMeshVert, saddleNeighborEdges, gridSize,
                                                                       contourEdges, edgesToRemove)

            saddleAllContourEdges.append(newEdges)
            saddleAllContourWeights.append(newWeights)
            saddleAllContourHeights.append(newHeights)

            # plot edges up to a certain step number
            # use this to check that the edges are in the correct order and not randomly jumping around
            # plotEdges(newEdges, 50)

        print("Num countour lines for saddle ", iNode, ":", len(saddleAllContourHeights))
        plotSaddleCountourLine(saddleAllContourEdges, saddleAllContourWeights, gridSize)
        # plt.show()
        plt.waitforbuttonpress()

    # currentEdge = [flatten2DIndex(20, 10, gridSize), flatten2DIndex(20, 11, gridSize)]
    # previousEdge = [flatten2DIndex(18, 11, gridSize), flatten2DIndex(20000, 10000, gridSize)]
    # print(isEdge(currentEdge, gridSize))
    # print(isEdge(previousEdge, gridSize))


    # square = findSquareFromEdge(previousEdge, currentEdge, gridSize)
    # print(square)


