import glob
import csv
import pandas as pd

from M01_TopologicalExtraction import *
from M03_Preprocessing import *

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
    ]

    inputSegmentedField = './Data/geodesics/JulienExample/tree1segs.vtk'
    inputMergeTreeNodes = './Data/geodesics/JulienExample/tree1nodes.vtk'
    inputMergeTreeEdges = './Data/geodesics/JulienExample/tree1edges.vtk'

    inputSegmentedField2 = './Data/geodesics/JulienExample/tree2segs.vtk'
    inputMergeTreeNodes2 = './Data/geodesics/JulienExample/tree2nodes.vtk'
    inputMergeTreeEdges2 = './Data/geodesics/JulienExample/tree2edges.vtk'

    # the 257 in main
    gridSize = (256, 257) #(Y, X) resolution in paraview (rows + 1, cols + 1)

    tree = Tree()
    tree.load(inputMergeTreeNodes, inputMergeTreeEdges, inputSegmentedField, gridSize=gridSize, splitTree=True, segmentationDataScalarName="multiGaussian0")
    tree.saddleTypeId = 1

    tree.extractContourLineConstraints()
    tree.reOrderContourline(True)


    for intermediateTreeDataFile in deformationDataFileToTree2s:
        intermediateData = pd.read_csv(intermediateTreeDataFile)

        preprocessTreeNodeData(tree, intermediateData)

        # rows = []
        # csvreader = csv.reader(file)
        # header = next(csvreader)
        # for row in csvreader:
        #     rows.append(row)
        # print(header)
        # print(rows)


    tree2 = Tree()
    tree2.load(inputMergeTreeNodes2, inputMergeTreeEdges2, inputSegmentedField2, gridSize=gridSize, splitTree=True, segmentationDataScalarName="multiGaussian1")
    tree2.saddleTypeId = 1

    tree2.extractContourLineConstraints()
    tree2.reOrderContourline(True)



