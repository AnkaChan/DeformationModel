import glob
import csv
import pandas as pd

from M01_TopologicalExtraction import *
from M03_Preprocessing import *
from M04_InterpolationAnimation import *
if __name__ == '__main__':
    # inDeformDataFolder = r'.\Data\geodesics\JulienExample\animation'
    inDeformDataFolder = r'.\Data\geodesics\JulienExample\pointData'

    nameToTree2s = "tree1ToInterpolated"
    deformationDataFileToTree2s = glob.glob(join(inDeformDataFolder, nameToTree2s + "*.csv"))

    treeToTreeCorrespondence = [
        [3, -1],
        [4, 2],
        [0, 0],
        [1, -1],
        [2, 1],
        [5, 3]
    ]

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

    animation = LinearAnimation(tree0, tree1, correspondences=treeToTreeCorrespondence)

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






