from S03_SolveScalarFieldContourLIneConstraintInterpolationConstraint import obj2vtkFolder
import pickle
import os


def saveList(outfile, listToSave):
    with open(outfile, 'wb') as fp:
        pickle.dump(listToSave, fp)

def readList(outfile):
    with open(outfile, 'rb') as fp:
        itemList = pickle.load(fp)
    return itemList

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


