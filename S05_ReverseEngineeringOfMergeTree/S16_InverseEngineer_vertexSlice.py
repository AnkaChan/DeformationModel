import json

from M02_ScalarFieldInterpolation import *
from M03_Preprocessing import *
from M05_InverseConstraints import *

if __name__ == '__main__':
    # inputSegmentedField = r'./Data/geodesics/VortexSlice/monoMesh_26_28/st/tree1segs.vtk'
    # inputMergeTreeNodes = r'./Data/geodesics/VortexSlice/monoMesh_26_28/st/tree1nodes.vtk'
    # inputMergeTreeEdges = r'./Data/geodesics/VortexSlice/monoMesh_26_28/st/tree1edges.vtk'
    #
    # inputSegmentedField2 = r'./Data/geodesics/VortexSlice/monoMesh_26_28/st/tree2segs.vtk'
    # inputMergeTreeNodes2 = r'./Data/geodesics/VortexSlice/monoMesh_26_28/st/tree2nodes.vtk'
    # inputMergeTreeEdges2 = r'./Data/geodesics/VortexSlice/monoMesh_26_28/st/tree2edges.vtk'

    inputSegmentedField = r'./Data/geodesics/HeatedFlowY/data_1_3/st/tree1segs.vtk'
    inputMergeTreeNodes = r'./Data/geodesics/HeatedFlowY/data_1_3/st/tree1nodes.vtk'
    inputMergeTreeEdges = r'./Data/geodesics/HeatedFlowY/data_1_3/st/tree1edges.vtk'

    inputSegmentedField2 = r'./Data/geodesics/HeatedFlowY/data_1_3/st/tree2segs.vtk'
    inputMergeTreeNodes2 = r'./Data/geodesics/HeatedFlowY/data_1_3/st/tree2nodes.vtk'
    inputMergeTreeEdges2 = r'./Data/geodesics/HeatedFlowY/data_1_3/st/tree2edges.vtk'

    # the 257 in main
    gridSize = (450, 150)
    # gridSize = (84, 212)
    # gridSize = (64, 192)  # (Y, X) resolution in paraview (rows + 1, cols + 1)
    # gridSize = (192, 64)  # (Y, X) resolution in paraview (rows + 1, cols + 1)

    tree0 = Tree()
    tree0.load(inputMergeTreeNodes, inputMergeTreeEdges, inputSegmentedField, gridSize=gridSize, splitTree=True,
               segmentationDataScalarName="velocityMagnitude")
    tree0.saddleTypeId = 2

    print("tree 0")
    tree0.extractContourLineConstraints3(draw=False)
    # tree0.extractContourLineConstraints()
    # tree0.reOrderContourline(False, waitTime=100)

    tree1 = Tree()
    tree1.load(inputMergeTreeNodes2, inputMergeTreeEdges2, inputSegmentedField2, gridSize=gridSize, splitTree=True,
               segmentationDataScalarName="velocityMagnitude")
    tree1.saddleTypeId = 2

    print("tree 1")
    tree1.extractContourLineConstraints3(draw=False)
    # tree1.reOrderContourline(False, waitTime=100)
