import json

from M02_ScalarFieldInterpolation import *
from M03_Preprocessing import *
from M04_InterpolationAnimation import *
from M05_InverseConstraints import *



if __name__ == '__main__':
    inputSegmentedField = r'.\Data\geodesics\MovingGaussian\monoMesh_2_3\st\tree1segs.vtk'
    inputMergeTreeNodes = r'.\Data\geodesics\MovingGaussian\monoMesh_2_3\st\tree1nodes.vtk'
    inputMergeTreeEdges = r'.\Data\geodesics\MovingGaussian\monoMesh_2_3\st\tree1edges.vtk'

    inputSegmentedField2 = r'.\Data\geodesics\MovingGaussian\monoMesh_2_3\st/tree2segs.vtk'
    inputMergeTreeNodes2 = r'.\Data\geodesics\MovingGaussian\monoMesh_2_3\st/tree2nodes.vtk'
    inputMergeTreeEdges2 = r'.\Data\geodesics\MovingGaussian\monoMesh_2_3\st/tree2edges.vtk'

    intermediateTreeFolder = r'./Data/geodesics/MovingGaussian/monoMesh_2_3/intermediateTree'

    intermediateTreeFiles = sorted(glob.glob(join(intermediateTreeFolder, "intermediateTree_*.csv")))
    intermediateTreeEdgesFiles = sorted(glob.glob(join(intermediateTreeFolder, "intermediateTreeEdge_*.csv")))


    # the 257 in main
    gridSize = (150, 150)  # (Y, X) resolution in paraview (rows + 1, cols + 1)
    # gridSize = (192, 64)  # (Y, X) resolution in paraview (rows + 1, cols + 1)

    tree0 = Tree()
    tree0.load(inputMergeTreeNodes, inputMergeTreeEdges, inputSegmentedField, gridSize=gridSize, splitTree=True,
               segmentationDataScalarName="Scalars_")
    tree0.saddleTypeId = 2

    tree0.extractContourLineConstraints3(draw=False)
    # tree0.extractContourLineConstraints()
    # tree0.reOrderContourline(True, waitTime=100)

    tree1 = Tree()
    tree1.load(inputMergeTreeNodes2, inputMergeTreeEdges2, inputSegmentedField2, gridSize=gridSize, splitTree=True,
               segmentationDataScalarName="Scalars_")
    tree1.saddleTypeId = 2
    tree1.extractContourLineConstraints3(draw=False)

    # tree1.extractContourLineConstraintsNew()
    # tree1.reOrderContourline(False, waitTime=100)

    treeMatcher = IntermediateTreeMatcher(tree0, tree1, intermediateTreeFiles, intermediateTreeEdgesFiles)
    
    

