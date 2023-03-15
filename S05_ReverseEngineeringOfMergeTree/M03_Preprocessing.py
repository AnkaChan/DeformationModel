from S03_SolveScalarFieldContourLIneConstraintInterpolationConstraint import obj2vtkFolder
import pickle
import os
import numpy as np

def saveList(outfile, listToSave):
    with open(outfile, 'wb') as fp:
        pickle.dump(listToSave, fp)

def readList(outfile):
    with open(outfile, 'rb') as fp:
        itemList = pickle.load(fp)
    return itemList

def preprocessTreeNodeData(tree1, intermediateTreeData):
    tree1Nodes = []
    intermediateTreeNodes = []

    # match Tree1 To nodes
    tree1NodesRef = []
    for node in tree1.nodes:
        tree1NodesRef.append(node.scalar)

    tree1NodesRef = np.array(tree1NodesRef)
    print(tree1NodesRef)
    print(intermediateTreeData['Points:1'].to_numpy())

    intermediateTreeDataHeights = intermediateTreeData['Points:1'].to_numpy()

    intermediateTreeDataHeights = (intermediateTreeDataHeights-np.min(intermediateTreeDataHeights)) * (np.max(tree1NodesRef) - np.min(tree1NodesRef)) / (np.max(intermediateTreeDataHeights) - np.min(intermediateTreeDataHeights))




if __name__ == "__main__":
    # generate the .vtk file of the terrain, adding the "height" column
    # obj2vtkFolder('./Data/S11_Geodesic.py/larger_domain/', outVtkFolder='./Data/S11_Geodesic.py/larger_domain/')
    # obj2vtkFolder('./Data/S11_Geodesic.py/larger_domain/MultiGaussian0/results/', outVtkFolder='./Data/S11_Geodesic.py/larger_domain/MultiGaussian0/results/')
    # obj2vtkFolder('./Data/S11_Geodesic.py/original/MultiGaussian0/', outVtkFolder='./Data/S11_Geodesic.py/original/MultiGaussian0/')
    # obj2vtkFolder('./Data/S11_Geodesic.py/original/MultiGaussian1/', outVtkFolder='./Data/S11_Geodesic.py/original/MultiGaussian1/')
    # obj2vtkFolder('./Data/S11_Geodesic.py/larger_domain/MultiGaussian0/', outVtkFolder='./Data/S11_Geodesic.py/larger_domain/MultiGaussian0/')
    contourInfoPath = './Data/S11_Geodesic.py/larger_domain/MultiGaussian0/ContourInfo'
    contourEdges = readList(os.path.join(contourInfoPath, "contourEdges.txt"))
    print(contourEdges)
    print(type(contourEdges))


