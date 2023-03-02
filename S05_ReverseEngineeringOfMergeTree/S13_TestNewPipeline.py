from M01_TopologicalExtraction import *

if __name__ == '__main__':
    # inputScalarField = './Data/S04_GenerateNewScalarField.py/M31562.obj'
    # inputSegmentedField = './Data/S04_GenerateNewScalarField.py/Topology/M31562/Seg.vtk'
    # inputMergeTreeNodes = './Data/S04_GenerateNewScalarField.py/Topology/M31562/Node.vtk'
    # inputMergeTreeEdges = './Data/S04_GenerateNewScalarField.py/Topology/M31562/Edges.vtk'
    # gridSize = (200, 200)

    # inputSegmentedField = './Data/WassasiteinDistance/tree1segs.vtk'
    # inputMergeTreeNodes = './Data/WassasiteinDistance/tree1nodes.vtk'
    # inputMergeTreeEdges = './Data/WassasiteinDistance/tree1edges.vtk'

    inputSegmentedField = './Data/geodesics/JulienExample/tree1segs.vtk'
    inputMergeTreeNodes = './Data/geodesics/JulienExample/tree1nodes.vtk'
    inputMergeTreeEdges = './Data/geodesics/JulienExample/tree1edges.vtk'


    # the 257 in main
    gridSize = (256, 257) #(Y, X) resolution in paraview (rows + 1, cols + 1)

    tree = Tree()
    tree.load(inputMergeTreeNodes, inputMergeTreeEdges, inputSegmentedField, gridSize=gridSize, splitTree=True, segmentationDataScalarName="multiGaussian0")
    tree.saddleTypeId = 1

    tree.extractContourLineConstraints()

    tree.reOrderContourline()

    # for iSaddle in range()
    # print("Num countour lines for saddle ", iNode, ":", len(saddleAllContourHeights))
    # plotSaddleCountourLine(saddleAllContourEdges, saddleAllContourWeights, gridSize)