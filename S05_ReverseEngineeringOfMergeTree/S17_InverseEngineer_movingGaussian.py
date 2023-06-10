import json

from M02_ScalarFieldInterpolation import *
from M03_Preprocessing import *
from M05_InverseConstraints import *

if __name__ == '__main__':
    inputSegmentedField = r'./Data/geodesics/MovingGaussian/monoMesh_0_2/mt/tree1segs.vtk'
    inputMergeTreeNodes = r'./Data/geodesics/MovingGaussian/monoMesh_0_2/mt/tree1nodes.vtk'
    inputMergeTreeEdges = r'./Data/geodesics/MovingGaussian/monoMesh_0_2/mt/tree1edges.vtk'

    inputSegmentedField2 = r'/Data/geodesics/MovingGaussian/monoMesh_0_2/mt/tree2segs.vtk'
    inputMergeTreeNodes2 = r'/Data/geodesics/MovingGaussian/monoMesh_0_2/mt/tree2nodes.vtk'
    inputMergeTreeEdges2 = r'/Data/geodesics/MovingGaussian/monoMesh_0_2/mt/tree2edges.vtk'

    # the 257 in main
    gridSize = (150, 150)  # (Y, X) resolution in paraview (rows + 1, cols + 1)
    # gridSize = (192, 64)  # (Y, X) resolution in paraview (rows + 1, cols + 1)

    tree0 = Tree()
    tree0.load(inputMergeTreeNodes, inputMergeTreeEdges, inputSegmentedField, gridSize=gridSize, splitTree=True,
               segmentationDataScalarName="Scalars_")
    tree0.saddleTypeId = 1

    tree0.extractContourLineConstraints3(draw=True)
    # tree0.extractContourLineConstraints()
    # tree0.reOrderContourline(True, waitTime=100)

    tree1 = Tree()
    tree1.load(inputMergeTreeNodes2, inputMergeTreeEdges2, inputSegmentedField2, gridSize=gridSize, splitTree=True,
               segmentationDataScalarName="Scalars_")
    tree1.saddleTypeId = 1

    tree1.extractContourLineConstraintsNew()
    tree1.reOrderContourline(False, waitTime=100)
