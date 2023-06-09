import json

from M02_ScalarFieldInterpolation import *
from M03_Preprocessing import *
from M05_InverseConstraints import *

if __name__ == '__main__':
    inputSegmentedField = r'./Data/geodesics/VortexSlice/monoMesh_8_10_smoothed/mt/tree1segs.vtk'
    inputMergeTreeNodes = r'./Data/geodesics/VortexSlice/monoMesh_8_10_smoothed/mt/tree1nodes.vtk'
    inputMergeTreeEdges = r'./Data/geodesics/VortexSlice/monoMesh_8_10_smoothed/mt/tree1edges.vtk'

    inputSegmentedField2 = r'./Data/geodesics/VortexSlice/monoMesh_8_10_smoothed/mt/tree2segs.vtk'
    inputMergeTreeNodes2 = r'./Data/geodesics/VortexSlice/monoMesh_8_10_smoothed/mt/tree2nodes.vtk'
    inputMergeTreeEdges2 = r'./Data/geodesics/VortexSlice/monoMesh_8_10_smoothed/mt/tree2edges.vtk'

    # the 257 in main
    gridSize = (64, 192)  # (Y, X) resolution in paraview (rows + 1, cols + 1)
    # gridSize = (192, 64)  # (Y, X) resolution in paraview (rows + 1, cols + 1)

    tree0 = Tree()
    tree0.load(inputMergeTreeNodes, inputMergeTreeEdges, inputSegmentedField, gridSize=gridSize, splitTree=True,
               segmentationDataScalarName="speed_smoothed")
    tree0.saddleTypeId = 1

    tree0.extractContourLineConstraints3()
    # tree0.extractContourLineConstraints()
    tree0.reOrderContourline(False, waitTime=100)

    tree1 = Tree()
    tree1.load(inputMergeTreeNodes2, inputMergeTreeEdges2, inputSegmentedField2, gridSize=gridSize, splitTree=True,
               segmentationDataScalarName="speed_smoothed")
    tree1.saddleTypeId = 1

    tree1.extractContourLineConstraintsNew()
    tree1.reOrderContourline(False, waitTime=100)
