from M01_TopologicalExtraction import *
from M02_ScalarFieldInterpolation import *

if __name__ == '__main__':
    # inSF = './Data/geodesics/VortexSlice/monoMesh_sub_8_10_Padded/monoMesh_008_smoothed.npy'
    # outSF = './Data/geodesics/VortexSlice/monoMesh_sub_8_10_Padded/monoMesh_008_smoothed_padded.npy'

    inSF = './Data/geodesics/VortexSlice/monoMesh_sub_8_10_Padded/monoMesh_010_smoothed.npy'
    outSF = './Data/geodesics/VortexSlice/monoMesh_sub_8_10_Padded/monoMesh_010_smoothed_padded.npy'

    paddingSize = 10

    sf = np.load(inSF)

    paddingBoundaryVal = np.min(sf)

    print(sf.shape)

    newShape = (sf.shape[0] + 2 * paddingSize, sf.shape[1] + 2 * paddingSize)

    newSf = np.zeros(newShape, dtype=sf.dtype)

    laplacianMatFile = "LaplacianMat_{:d}_{:d}.npz".format(newShape[0], newShape[1])

    dType = np.float64

    if os.path.exists(laplacianMatFile):
        P = sparse.load_npz(laplacianMatFile).astype(dType)
    else:
        P = generateLaplacianMat(newShape)
        sparse.save_npz(laplacianMatFile, P)


    xRange = [-2, 2]
    yRange = [-2, 2]
    X = np.linspace(*xRange, newShape[1])
    Y = np.linspace(*yRange, newShape[0])
    X, Y = np.meshgrid(X, Y)

    # boundary:


    weights = []
    verts = []
    val = []
    # for i in range(newShape[0]):
    #     weights.append((1.0,))
    #     verts.append((flatten2DIndex(i, 0, newShape),))
    #     val.append(paddingBoundaryVal)

        # weights.append((1.0,))
        # verts.append((flatten2DIndex(i, newShape[1]-1, newShape),))
        # val.append(paddingBoundaryVal)

    # for i in range(1, newShape[1]):
    #     weights.append((1.0,))
    #     verts.append((flatten2DIndex(0, i, newShape),))
    #     val.append(paddingBoundaryVal)
    #
    #     weights.append((1.0,))
    #     verts.append((flatten2DIndex(newShape[0] - 1, i, newShape),))
    #     val.append(paddingBoundaryVal)

    for i in range(sf.shape[0]):
        for j in range(sf.shape[1]):
            weights.append((1.0,))
            verts.append((flatten2DIndex(i+paddingSize, j+paddingSize, newShape),))
            val.append(sf[i,j])

    field = interpolateScalarField(newShape, P, verts,
                                   weights, val,
                                   [], [], [],
                                   fixBoundary=True, boundaryValue=paddingBoundaryVal,
                                   dType=np.float64, )

    padded = field.reshape(*newShape)

    assert np.max(np.abs(padded[paddingSize:-paddingSize, paddingSize:-paddingSize]-sf)) < 1e-4

    np.save(outSF, padded, )

    outFile = join("TestPadding.obj")
    writeOBj(outFile, X, Y, field, newShape)