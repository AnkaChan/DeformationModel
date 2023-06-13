from M01_TopologicalExtraction import *
from M02_ScalarFieldInterpolation import *

import pyvista as pv
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from pyevtk.hl import imageToVTK

if __name__ == '__main__':
    # inSF = './Data/geodesics/VortexSlice/monoMesh_sub_8_10_Padded/monoMesh_008_smoothed.npy'
    # outSF = './Data/geodesics/VortexSlice/monoMesh_sub_8_10_Padded/monoMesh_008_smoothed_padded.npy'

    # inSF = './Data/geodesics/VortexSlice/monoMesh_sub_8_10_Padded/monoMesh_010_smoothed.npy'
    # outSF = './Data/geodesics/VortexSlice/monoMesh_sub_8_10_Padded/monoMesh_010_smoothed_padded.npy'

    input_dir = "Data/geodesics/VortexSlice/"
    # output_dir = "./VortexSlice/monoMesh_8_10/"
    inputFile = "monoMesh_030.vti"
    outfile = input_dir + "monoMesh_030_smoothed"
    scalar_name = 'speed'
    new_scalar_name = 'speed_smoothed'
    gridSize = (192, 64)
    smoothing = True
    flipping = True
    padding = True

    inputSf = pv.read(input_dir + inputFile)
    sf = np.array(inputSf[scalar_name]).reshape((gridSize[1], gridSize[0])).T

    if smoothing:
        sf = gaussian_filter(sf, sigma=1)
        smoothed_pts = sf.reshape(gridSize[0], gridSize[1], 1)

    if flipping:
        sf = -1 * sf
        smoothed_pts = sf.reshape(gridSize[0], gridSize[1], 1)

    if padding:
        paddingSize = 10
        # sf = np.load(inSF)
        paddingBoundaryVal = np.min(sf)
        print(sf.shape)
        newShape = (sf.shape[0] + 2 * paddingSize, sf.shape[1] + 2 * paddingSize)
        # newSf = np.zeros(newShape, dtype=sf.dtype)

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
        smoothed_pts = padded.reshape(newShape[0], newShape[1], 1)
        print(padded.shape)
        plt.imshow(padded)
        plt.show()
        assert np.max(np.abs(padded[paddingSize:-paddingSize, paddingSize:-paddingSize]-sf)) < 1e-4

    print(outfile)
    imageToVTK(outfile, pointData={"speed_smoothed": smoothed_pts})
    # np.save(outSF, padded, )
    #
    # outFile = join("TestPadding.obj")
    # writeOBj(outFile, X, Y, field, newShape)