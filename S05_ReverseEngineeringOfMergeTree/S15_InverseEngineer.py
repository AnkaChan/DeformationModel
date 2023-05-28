import json

from M02_ScalarFieldInterpolation import *
from M03_Preprocessing import *
from M05_InverseConstraints import *

if __name__ == '__main__':

    outFolder = r'./Data/' + os.path.basename(__file__)
    outFolderConstraints = r'./Data/' + os.path.basename(__file__) + "/Constraints"

    os.makedirs(outFolder, exist_ok=True)
    os.makedirs(outFolderConstraints, exist_ok=True)

    inputFolder = r'Data\S14_TreeDeformation_1.py'

    allFiles = glob.glob(join(inputFolder, "*.json"))

    constraintData0 = json.load(open(allFiles[0]))
    gridSize = constraintData0["gridSize"]
    laplacianMatFile = "LaplacianMat_{:d}_{:d}.npz".format(gridSize[0], gridSize[1])

    dType = np.float64

    if os.path.exists(laplacianMatFile):
        P = sparse.load_npz(laplacianMatFile).astype(dType)
    else:
        P = generateLaplacianMat(gridSize)
        sparse.save_npz(laplacianMatFile, P)

    intersector = ContourLineGridIntersector(gridSize)

    xRange = [-2, 2]
    yRange = [-2, 2]
    X = np.linspace(*xRange, gridSize[1])
    Y = np.linspace(*yRange, gridSize[0])
    X, Y = np.meshgrid(X, Y)
    # writeOBj(outFile, 200, X, Y, Z)
    for inFile in tqdm.tqdm(allFiles):
        fileName = Path(inFile).stem
        constraintData = json.load(open(inFile))

        intersector.intersectWithGrid(constraintData, plot=True)
        constraintFile = join(outFolderConstraints, fileName + '.ply')
        intersector.saveConstraints(constraintFile, xRange, yRange)

        field = interpolateScalarField(gridSize, P, intersector.equalityConstraintIds, intersector.equalityConstraintWeights, intersector.equalityConstraintVal,
                               [], [], [],
                               fixBoundary=True, boundaryValue=0,
                               dType=np.float64,)

        outFile = join(outFolder, fileName + ".obj")
        writeOBj(outFile, X, Y, field, gridSize)

    obj2vtkFolder(outFolder,)
