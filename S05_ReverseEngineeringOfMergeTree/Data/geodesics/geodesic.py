#!/usr/bin/env python

# /// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
# /// \date February 2022.

# uncomment the following three lines to ensure this script works in future versions
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 10

from paraview.simple import *

# ------------------------------------------------------------------------------
# 1) Creating two Gaussian mixtures, multiGaussian0 and multiGaussian1
# ------------------------------------------------------------------------------

# create a new 'Plane'
plane1 = Plane(registrationName='Plane1')
plane1.XResolution = 256
plane1.YResolution = 255
plane1.Origin = [-0.7, -0.7, 0.0]
plane1.Point1 = [0.7, -0.7, 0.0]
plane1.Point2 = [-0.7, 0.7, 0.0]

# create a new 'Python Calculator'
pythonCalculator1 = PythonCalculator(registrationName='PythonCalculator1', Input=plane1)
pythonCalculator1.Expression = '-exp(-((inputs[0].Points[:,0]-0.25)*(inputs[0].Points[:,0]-0.25) + (inputs[0].Points[:,1]-0.25)*(inputs[0].Points[:,1]-0.25))/(2*0.2*0.2))'
pythonCalculator1.ArrayName = 'gaussian0'

# create a new 'Python Calculator'
pythonCalculator2 = PythonCalculator(registrationName='PythonCalculator2', Input=pythonCalculator1)
pythonCalculator2.Expression = '-exp(-((inputs[0].Points[:,0]+0.25)*(inputs[0].Points[:,0]+0.25) + (inputs[0].Points[:,1]+0.25)*(inputs[0].Points[:,1]+0.25))/(2*0.2*0.2))'
pythonCalculator2.ArrayName = 'gaussian1'

# create a new 'Python Calculator'
pythonCalculator3 = PythonCalculator(registrationName='PythonCalculator3', Input=pythonCalculator2)
pythonCalculator3.Expression = '-exp(-((inputs[0].Points[:,0]+0.15)*(inputs[0].Points[:,0]+0.15) + (inputs[0].Points[:,1]-0.25)*(inputs[0].Points[:,1]-0.25))/(2*0.1*0.1))'
pythonCalculator3.ArrayName = 'gaussian2'

# create a new 'Python Calculator'
pythonCalculator4 = PythonCalculator(registrationName='PythonCalculator4', Input=pythonCalculator3)
pythonCalculator4.Expression = '1.5*gaussian0+1.2*gaussian1+0.85*gaussian2'
pythonCalculator4.ArrayName = 'multiGaussian0'

# create a new 'Python Calculator'
pythonCalculator5 = PythonCalculator(registrationName='PythonCalculator5', Input=pythonCalculator4)
pythonCalculator5.Expression = '1.5*gaussian0+1.2*gaussian1'
pythonCalculator5.ArrayName = 'multiGaussian1'

# create a new 'Tetrahedralize'
tetrahedralize1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=pythonCalculator5)

# ------------------------------------------------------------------------------
# 2) Computing the merge trees
# ------------------------------------------------------------------------------
output_dir = "./Data/geodesics/JulienExample/"

# create a new 'TTK Merge and Contour Tree (FTM)'
tTKMergeandContourTreeFTM1 = TTKMergeandContourTreeFTM(registrationName='TTKMergeandContourTreeFTM1', Input=tetrahedralize1)
tTKMergeandContourTreeFTM1.ScalarField = ['POINTS', 'multiGaussian0']
tTKMergeandContourTreeFTM1.InputOffsetField = ['POINTS', 'gaussian0']
tTKMergeandContourTreeFTM1.TreeType = 'Join Tree'
# save tree 1 info
SetActiveSource(tTKMergeandContourTreeFTM1)
SaveData(output_dir + 'tree1nodes.vtk', proxy=tTKMergeandContourTreeFTM1, PointDataArrays=['CriticalType', 'NodeId', 'RegionSize', 'RegionSpan', 'Scalar', 'VertexId'])
SetActiveSource(tTKMergeandContourTreeFTM1)
tTKMergeandContourTreeFTM1_1 = GetActiveSource()
SaveData(output_dir + 'tree1edges.vtk', proxy=OutputPort(tTKMergeandContourTreeFTM1_1, 1), PointDataArrays=['Scalar', 'ttkMaskScalarField'],
    CellDataArrays=['RegionSize', 'RegionSpan', 'SegmentationId', 'downNodeId', 'upNodeId'])
SetActiveSource(tTKMergeandContourTreeFTM1)
tTKMergeandContourTreeFTM1_2 = GetActiveSource()
SaveData(output_dir + 'tree1segs.vtk', proxy=OutputPort(tTKMergeandContourTreeFTM1_2, 2), PointDataArrays=['Normals', 'RegionSize', 'RegionSpan', 'RegionType', 'SegmentationId', 'TextureCoordinates', 'gaussian0', 'gaussian1', 'gaussian2', 'multiGaussian0', 'multiGaussian1'])


# create a new 'TTK Merge and Contour Tree (FTM)'
tTKMergeandContourTreeFTM2 = TTKMergeandContourTreeFTM(registrationName='TTKMergeandContourTreeFTM2', Input=tetrahedralize1)
tTKMergeandContourTreeFTM2.ScalarField = ['POINTS', 'multiGaussian1']
tTKMergeandContourTreeFTM2.InputOffsetField = ['POINTS', 'gaussian0']
tTKMergeandContourTreeFTM2.TreeType = 'Join Tree'
# save tree 2 info
tTKMergeandContourTreeFTM2 = FindSource('TTKMergeandContourTreeFTM2')
SetActiveSource(tTKMergeandContourTreeFTM2)
SaveData(output_dir + 'tree2nodes.vtk', proxy=tTKMergeandContourTreeFTM2, PointDataArrays=['CriticalType', 'NodeId', 'RegionSize', 'RegionSpan', 'Scalar', 'VertexId'])
SetActiveSource(tTKMergeandContourTreeFTM2)
tTKMergeandContourTreeFTM2_1 = GetActiveSource()
SaveData(output_dir + 'tree2edges.vtk', proxy=OutputPort(tTKMergeandContourTreeFTM2_1, 1), PointDataArrays=['Scalar', 'ttkMaskScalarField'],
    CellDataArrays=['RegionSize', 'RegionSpan', 'SegmentationId', 'downNodeId', 'upNodeId'])
SetActiveSource(tTKMergeandContourTreeFTM2)
tTKMergeandContourTreeFTM2_2 = GetActiveSource()
SaveData(output_dir + 'tree2segs.vtk', proxy=OutputPort(tTKMergeandContourTreeFTM2_2, 2), PointDataArrays=['Normals', 'RegionSize', 'RegionSpan', 'RegionType', 'SegmentationId', 'TextureCoordinates', 'gaussian0', 'gaussian1', 'gaussian2', 'multiGaussian0', 'multiGaussian1'])


# ------------------------------------------------------------------------------
# 3) Grouping together the merge tree outputs for the Wasserstein distance
# ------------------------------------------------------------------------------

# find source
tTKMergeandContourTreeFTM2_1 = FindSource('TTKMergeandContourTreeFTM2')

# create a new 'Group Datasets'
groupDatasets2 = GroupDatasets(registrationName='GroupDatasets2', Input=[tTKMergeandContourTreeFTM2, OutputPort(tTKMergeandContourTreeFTM2_1,1)])
groupDatasets2.BlockNames = ['TTKMergeandContourTreeFTM2', 'TTKMergeandContourTreeFTM2']

# find source
tTKMergeandContourTreeFTM1_1 = FindSource('TTKMergeandContourTreeFTM1')

# create a new 'Group Datasets'
groupDatasets1 = GroupDatasets(registrationName='GroupDatasets1', Input=[tTKMergeandContourTreeFTM1, OutputPort(tTKMergeandContourTreeFTM1_1,1)])
groupDatasets1.BlockNames = ['TTKMergeandContourTreeFTM1', 'TTKMergeandContourTreeFTM1']

# create a new 'Group Datasets'
groupDatasets3 = GroupDatasets(registrationName='GroupDatasets3', Input=[groupDatasets1, groupDatasets2])
# group of merge trees -- input
groupDatasets3.BlockNames = ['GroupDatasets1', 'GroupDatasets2']

# ------------------------------------------------------------------------------
# 4) Computing the Wasserstein distance
# ------------------------------------------------------------------------------

# create a new 'TTK MergeTreeClustering'
# tTKMergeTreeClustering1 = TTKMergeTreeClustering(registrationName='TTKMergeTreeClustering1', Input=groupDatasets3,
#     OptionalInputclustering=None)
# tTKMergeTreeClustering1.ComputeBarycenter = 1
# # change this parameter to change the position of the interpolated tree
# # change to 1 to retrieve the final position of all nodes (including destroyed
# # nodes)
# tTKMergeTreeClustering1.Alpha = 0.44
# tTKMergeTreeClustering1.DimensionToshift = 'Z'
# tTKMergeTreeClustering1.Barycenterpositionaccordingtoalpha = 1
# tTKMergeTreeClustering1.ImportantPairs = 0.0
# tTKMergeTreeClustering1.ImportantPairsSpacing = 0.25

a = list(range(1, 100))
t_list = [i/100 for i in a]

for t in t_list:
    tTKMergeTreeClustering1 = TTKMergeTreeClustering(registrationName='TTKMergeTreeClustering1', Input=groupDatasets3,
                                                     OptionalInputclustering=None)
    tTKMergeTreeClustering1.ComputeBarycenter = 1
    tTKMergeTreeClustering1.Alpha = t
    tTKMergeTreeClustering1.DimensionToshift = 'Z'
    tTKMergeTreeClustering1.Barycenterpositionaccordingtoalpha = 1
    tTKMergeTreeClustering1.ImportantPairs = 0.0
    tTKMergeTreeClustering1.ImportantPairsSpacing = 0.25

    # ------------------------------------------------------------------------------
    # 5) Isolating the different outputs
    # ------------------------------------------------------------------------------

    # find source
    tTKMergeTreeClustering1_2 = FindSource('TTKMergeTreeClustering1')

    # create a new 'Extract Block'
    extractBlock4 = ExtractBlock(registrationName='ExtractBlock4', Input=OutputPort(tTKMergeTreeClustering1_2,2))
    extractBlock4.Selectors = ['/Root/Block1']

    # create a new 'Extract Block'
    extractBlock3 = ExtractBlock(registrationName='ExtractBlock3', Input=OutputPort(tTKMergeTreeClustering1_2,2))
    extractBlock3.Selectors = ['/Root/Block0']

    SaveData(output_dir + "animation/tree1ToInterpolated_" + str(int(t*100)) + ".csv", extractBlock3)
    SaveData(output_dir + "animation/interpolatedTotree2_" + str(int(t*100)) + ".csv", extractBlock4)

# # create a new 'Extract Block'
# extractBlock2 = ExtractBlock(registrationName='ExtractBlock2', Input=OutputPort(tTKMergeTreeClustering1_2, 0))
# extractBlock2.Selectors = ['/Root/Block0/Block0']
#
# # create a new 'Extract Block'
# extractBlock1 = ExtractBlock(registrationName='ExtractBlock1', Input=OutputPort(tTKMergeTreeClustering1_2, 0))
# extractBlock1.Selectors = ['/Root/Block1/Block0']


# SaveData("tree1ToInterpolated.csv", extractBlock3)
# SaveData("interpolatedTotree2.csv", extractBlock4)
#
# SaveData("nodesTree1.csv", extractBlock2, FieldAssociation='Point Data')
# SaveData("nodesTree2.csv", extractBlock1, FieldAssociation='Point Data')
