### Inverse Merge Tree
02/20/2023

__Geodesic.py__
This is a simple example with multigaussians. The script reverse engineers the input scalar field.
Input exmple:
+ inputScalarField: Data/Geodesic.py/larger_domain/MultiGaussian1/MultiGaussian1.obj (any .obj in the Geodesic.py folder works)
+ inputSegmentFeidl: Data/Geodesic.py/larger_domain/MultiGaussian1/Seg1.vtk
+ inputMergeTreeNodesField: Data/Geodesic.py/larger_domain/MultiGaussian1/Node1.vtk

Output:
Output folder: Data/Geodesic.py/larger_domain/MultiGaussian1/results/
Files in the output folder:
+ contourLine.ply: the contour line constraint we generated, the reverse engineered scalar field should mostly satisfy this constraint
+ result.obj: the reverse engineered scalar field
+ ScalarFieldSimplified.vtk: simplifed version of the original scalar field. In this example, there's not much of a difference between the original and simplified

__test.py__
This is a more complicated example using randomly generated scalar field from the script "S04_GeneratedNewScalarField.py".
The input files are similar from the above Geodesic.py example, the difference is that the input folder.
Input example:
+ inputScalarField: Data/S04_GenerateNewScalarField.py/M31562.obj
+ inputSegmentedField: Data/S04_GenerateNewScalarField.py/Topology/M31562/Seg.vtk
+ inputMergeTreeNodesField: Data/S04_GenerateNewScalarField.py/Topology/M31562/Node.vtk

I increased the number of directions in the `directions` variable to fix corner cases of the `findContourLineHeight()` function in the `M01_TopologicalExtraction.py` file. I think the high-level idea is that sometime TTK generates multiple small segments around the saddle and we need to find the correct segments to determine the contour line height. I increased the number of directions and sorted the segments based on how many of these directions land in each saddle. I took the three segments that appeared the most, in order to determine the height of the contour. This part of the code can be automated.

__ParameterizeContourLine.py__
This script takes the contour line information, traces the contour line from the saddle and back to the saddle. The original contour edge information is in a random order.
Input includes all the scalar field and merge tree information, and contour line information.

Ouput:
+ newEdges: a list of reordered contour edges
+ newWeights: a list of reordered contour weights
+ newHeights: a list of reordered contour heights

We can use the `plotEdges()` function to check that the edges are in the correct order.
`plotEdges(contourEdges, step_number)` plots the edges end points one by one up to a step number. If we increment the step number one by one, we will see the progression of the contour edges in the correct order around the saddle.