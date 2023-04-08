import glob
import csv
import pandas as pd

from M01_TopologicalExtraction import *
from M03_Preprocessing import *
from M04_InterpolationAnimation import *


def getIntermediateTreeNodeHeight(tree0, tree1, treeToTreeCorrespondence, nSteps = 100):
    h0 = np.zeros(tree0.numNodes())

    heights = np.zeros((tree0.numNodes(), nSteps))

    for iNode in range(tree0.numNodes()):
        assert tree0.node(iNode).id == iNode
        h0[iNode] = tree0.node(iNode).scalar
        corrTree1 = findCorr(iNode, 0, treeToTreeCorrespondence)

        if corrTree1 != -1:
            heights[iNode, :] = np.linspace(h0[iNode], tree1.node(corrTree1).scalar, nSteps, endpoint=True)
        # else:


    # find the higher node and the lower node
    # node 3
    lowerNodeIdTree0 = tree0.node(3).downNodes[0]
    lowerNodeIdTree1 = findCorr(lowerNodeIdTree0, 0, treeToTreeCorrespondence)
    assert lowerNodeIdTree1 != -1

    heightlow = tree1.node(lowerNodeIdTree1).scalar

    heighthigh = -1
    vanishingUpNode = 1

    for upNodeId in tree0.node(3).upNodes:
        upNodeIdCorrTree1 = findCorr(lowerNodeIdTree0, 0, treeToTreeCorrespondence)

        if findCorr(lowerNodeIdTree0, 0, treeToTreeCorrespondence) != -1:
            heighthigh = findCorr(lowerNodeIdTree0, 0, treeToTreeCorrespondence)

        else:
            vanishingUpNode = upNodeId

    heights[3, :] = np.linspace(h0[3], (heighthigh + heightlow) * 0.5, nSteps, endpoint=True)
    heights[vanishingUpNode, :] = np.linspace(h0[vanishingUpNode], (heighthigh + heightlow) * 0.5, nSteps, endpoint=True)

    return heights

if __name__ == '__main__':
    # inDeformDataFolder = r'.\Data\geodesics\JulienExample\animation'
    inDeformDataFolder = r'.\Data\geodesics\JulienExample\pointData'

    nameToTree2s = "tree1ToInterpolated"
    deformationDataFileToTree2s = glob.glob(join(inDeformDataFolder, nameToTree2s + "*.csv"))

    treeToTreeCorrespondence = np.array([
        [3, -1],
        [4, 2],
        [0, 0],
        [1, -1],
        [2, 1],
        [5, 3]
    ], dtype=np.int)



    inputSegmentedField = './Data/geodesics/JulienExample/tree1segs.vtk'
    inputMergeTreeNodes = './Data/geodesics/JulienExample/tree1nodes.vtk'
    inputMergeTreeEdges = './Data/geodesics/JulienExample/tree1edges.vtk'

    inputSegmentedField2 = './Data/geodesics/JulienExample/tree2segs.vtk'
    inputMergeTreeNodes2 = './Data/geodesics/JulienExample/tree2nodes.vtk'
    inputMergeTreeEdges2 = './Data/geodesics/JulienExample/tree2edges.vtk'

    # the 257 in main
    gridSize = (256, 257) #(Y, X) resolution in paraview (rows + 1, cols + 1)

    tree0 = Tree()
    tree0.load(inputMergeTreeNodes, inputMergeTreeEdges, inputSegmentedField, gridSize=gridSize, splitTree=True, segmentationDataScalarName="multiGaussian0")
    tree0.saddleTypeId = 1

    tree0.extractContourLineConstraintsNew()
    # tree0.extractContourLineConstraints()
    tree0.reOrderContourline(False, waitTime=100)

    tree1 = Tree()
    tree1.load(inputMergeTreeNodes2, inputMergeTreeEdges2, inputSegmentedField2, gridSize=gridSize, splitTree=True, segmentationDataScalarName="multiGaussian1")
    tree1.saddleTypeId = 1

    tree1.extractContourLineConstraintsNew()
    tree1.reOrderContourline(False, waitTime=100)

    intermediateTreeNodeHeight = getIntermediateTreeNodeHeight(tree0, tree1, treeToTreeCorrespondence,  nSteps = 100)

    animation = LinearAnimation(tree0, tree1, correspondences=treeToTreeCorrespondence, intermediateTreeHeights=intermediateTreeNodeHeight)


    # for intermediateTreeDataFile in deformationDataFileToTree2s:
    #     intermediateData = pd.read_csv(intermediateTreeDataFile)
    #
    #     preprocessTreeNodeData(tree, intermediateData)
    #
    #     # rows = []
    #     # csvreader = csv.reader(file)
    #     # header = next(csvreader)
    #     # for row in csvreader:
    #     #     rows.append(row)
    #     # print(header)
    #     # print(rows)






