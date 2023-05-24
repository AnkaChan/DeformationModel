# uncomment the following three lines to ensure this script works in future versions
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 10

from paraview.simple import *
import sys
import os

t = float(sys.argv[1])
print(t)
input_dir = "./HeatedFlowY/"
output_dir = "./HeatedFlowY/data_604_606/"

inputStartFile = "data_604.vti"
inputEndFile = "data_606.vti"
threshold = 0.06
treeType = "mt"
scalarField = 'velocityMagnitude'

mt_dir = os.path.join(output_dir, treeType)
print(mt_dir)

if not os.path.exists(mt_dir):
    os.makedirs(mt_dir)

inputStartData = XMLImageDataReader(FileName=[input_dir + inputStartFile])
inputEndData = XMLImageDataReader(FileName=[input_dir + inputEndFile])

# ------------------------------------------------------------------------------
# 0. Perform persistence Simplification
# ------------------------------------------------------------------------------

simplifiedStart = TTKTopologicalSimplificationByPersistence(registrationName='simplifiedStart', Input=inputStartData)
simplifiedStart.InputArray = ['POINTS', scalarField]
simplifiedStart.PersistenceThreshold = threshold

simplifiedEnd = TTKTopologicalSimplificationByPersistence(registrationName='simplifiedEnd', Input=inputEndData)
simplifiedEnd.InputArray = ['POINTS', scalarField]
simplifiedEnd.PersistenceThreshold = threshold

# ------------------------------------------------------------------------------
# 1. Compute and save start and end merge trees
# ------------------------------------------------------------------------------

# compute start merge trees
MTStart = TTKMergeandContourTreeFTM(registrationName='MTStart', Input=simplifiedStart)
MTStart.ScalarField = ['POINTS', scalarField]
MTStart.TreeType = 'Join Tree'

# save start tree info
SetActiveSource(MTStart)
SaveData(mt_dir + '/' + 'tree1nodes.vtk', proxy=MTStart, PointDataArrays=['CriticalType', 'NodeId', 'RegionSize', 'RegionSpan', 'Scalar', 'VertexId'])
SetActiveSource(MTStart)
MTStart_1 = GetActiveSource()
SaveData(mt_dir + '/' + 'tree1edges.vtk', proxy=OutputPort(MTStart_1, 1), PointDataArrays=['Scalar', 'ttkMaskScalarField'],
    CellDataArrays=['RegionSize', 'RegionSpan', 'SegmentationId', 'downNodeId', 'upNodeId'])
SetActiveSource(MTStart)
MTStart_2 = GetActiveSource()
SaveData(mt_dir + '/' + 'tree1segs.vtk', proxy=OutputPort(MTStart_2, 2))

# compute end merge tree
MTEnd = TTKMergeandContourTreeFTM(registrationName='MTEnd', Input=simplifiedEnd)
MTEnd.ScalarField = ['POINTS', scalarField]
MTEnd.TreeType = 'Join Tree'

# save end tree info
SetActiveSource(MTEnd)
SaveData(mt_dir + '/' + 'tree2nodes.vtk', proxy=MTEnd, PointDataArrays=['CriticalType', 'NodeId', 'RegionSize', 'RegionSpan', 'Scalar', 'VertexId'])
SetActiveSource(MTEnd)
MTEnd_1 = GetActiveSource()
SaveData(mt_dir + '/' + 'tree2edges.vtk', proxy=OutputPort(MTEnd_1, 1), PointDataArrays=['Scalar', 'ttkMaskScalarField'],
    CellDataArrays=['RegionSize', 'RegionSpan', 'SegmentationId', 'downNodeId', 'upNodeId'])
SetActiveSource(MTEnd)
MTEnd_2 = GetActiveSource()
SaveData(mt_dir + '/' + 'tree2segs.vtk', proxy=OutputPort(MTEnd_2, 2))

# ------------------------------------------------------------------------------
# 2. Grouping together the merge tree outputs for the Wasserstein distance
# ------------------------------------------------------------------------------

# find source
MTEnd_1 = FindSource('MTEnd')

# create a new 'Group Datasets'
groupDatasets2 = GroupDatasets(registrationName='GroupDatasets2', Input=[MTEnd, OutputPort(MTEnd_1,1)])
groupDatasets2.BlockNames = ['MTEnd', 'MTEnd']

# find source
MTStart_1 = FindSource('MTStart')

# create a new 'Group Datasets'
groupDatasets1 = GroupDatasets(registrationName='GroupDatasets1', Input=[MTStart, OutputPort(MTStart_1,1)])
groupDatasets1.BlockNames = ['MTStart', 'MTStart']

# create a new 'Group Datasets'
MTGroup = GroupDatasets(registrationName='MTGroup', Input=[groupDatasets1, groupDatasets2])
# group of merge trees -- input
MTGroup.BlockNames = ['MTStartGroup', 'MTEndGroup']

# ------------------------------------------------------------------------------
# 3. Computing the Wasserstein distance
# ------------------------------------------------------------------------------

# create a new 'TTK MergeTreeClustering'
MTClustering = TTKMergeTreeClustering(registrationName='MTClustering', Input=MTGroup,
    OptionalInputclustering=None)
MTClustering.ComputeBarycenter = 1
# change this parameter to change the position of the interpolated tree
# change to 1 to retrieve the final position of all nodes (including destroyed
# nodes)
# t = 0.45
MTClustering.Alpha = t
MTClustering.DimensionToshift = 'Z'
MTClustering.Barycenterpositionaccordingtoalpha = 1
MTClustering.ImportantPairs = 0.0
MTClustering.ImportantPairsSpacing = 0.25

# ------------------------------------------------------------------------------
# 5) Isolating the different outputs
# ------------------------------------------------------------------------------

# find source
MTClustering_2 = FindSource('MTClustering')

# create a new 'Extract Block'
extractBlock4 = ExtractBlock(registrationName='ExtractBlock4', Input=OutputPort(MTClustering_2,2))
extractBlock4.Selectors = ['/Root/Block1']

# create a new 'Extract Block'
extractBlock3 = ExtractBlock(registrationName='ExtractBlock3', Input=OutputPort(MTClustering_2,2))
extractBlock3.Selectors = ['/Root/Block0']

matchings_dir = os.path.join(output_dir, "matchings")

if not os.path.exists(matchings_dir):
    os.makedirs(matchings_dir)

# save matchings
SaveData(matchings_dir + "/tree1ToInterpolated_" + str(int(t*100)) + ".csv", extractBlock3)
SaveData(matchings_dir + "/interpolatedTotree2_" + str(int(t*100)) + ".csv", extractBlock4)

intermediate_dir = os.path.join(output_dir, "intermediateTree")

if not os.path.exists(intermediate_dir):
    os.makedirs(intermediate_dir)

# get the intermediate tree data
MTClustering = FindSource('MTClustering')
SetActiveSource(MTClustering)
intermediateTree = GetActiveSource()
SaveData(intermediate_dir + '/intermediateTree_' + str(int(t*100)) + ".csv", proxy=OutputPort(intermediateTree, 1),
    PointDataArrays=['BranchNodeID', 'ClusterID', 'CriticalType', 'NodeId', 'PercentMatchNode', 'Persistence', 'Scalar', 'TreeID', 'VertexId', 'isDummyNode', 'isImportantPair'],
    CellDataArrays=['BranchID', 'ClusterID', 'PercentMatchArc', 'Persistence', 'TreeID', 'downNodeId', 'isDummyArc', 'isImportantPair', 'upNodeId'],
    FieldDataArrays=['ClusterAssignment'],
    AddMetaData=0)
SaveData(intermediate_dir + '/intermediateTreeEdge_' + str(int(t*100)) + ".csv", proxy=OutputPort(intermediateTree, 1),
    PointDataArrays=['BranchNodeID', 'ClusterID', 'CriticalType', 'NodeId', 'PercentMatchNode', 'Persistence', 'Scalar', 'TreeID', 'VertexId', 'isDummyNode', 'isImportantPair'],
    CellDataArrays=['BranchID', 'ClusterID', 'PercentMatchArc', 'Persistence', 'TreeID', 'downNodeId', 'isDummyArc', 'isImportantPair', 'upNodeId'],
    FieldDataArrays=['ClusterAssignment'],
    FieldAssociation='Cell Data')
