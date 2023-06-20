import pyvista as pv
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from pyevtk.hl import imageToVTK

from ..M01_TopologicalExtraction import *
from ..M02_ScalarFieldInterpolation import *

input_dir = "./VortexSlice/"
# output_dir = "./VortexSlice/monoMesh_8_10/"
inputStartFile = "monoMesh_008.vti"
outfile = input_dir + "monoMesh_008_smoothed"
gridSize = (192, 64)

# smooth with Gaussian filter
fig = plt.figure()
input = pv.read(input_dir + inputStartFile)
# input.plot()
input_image = np.array(input['speed']).reshape((gridSize[1], gridSize[0])).T
ax1 = fig.add_subplot(121)  # left side
ax2 = fig.add_subplot(122)  # right side
ax1.imshow(input_image)
smoothed = gaussian_filter(input_image, sigma=1)
print(smoothed)
ax2.imshow(smoothed)
plt.show()
print(input_image.shape)
# speed_pts = input_image.reshape(gridSize[0], gridSize[1], 1)
# smoothed_pts = smoothed.reshape(gridSize[0], gridSize[1], 1)
print(outfile)

# pad the scalar field so the contours on the boundaries are not broken
paddingSize = 10

newShape = (smoothed.shape[0] + 2 * paddingSize, smoothed.shape[1] + 2 * paddingSize)
paddedSf = np.zeros(newShape, dtype=smoothed.dtype)
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

weights = []
verts = []
val = []

for i in range(smoothed.shape[0]):
    for j in range(smoothed.shape[1]):
        weights.append((1.0,))
        verts.append((flatten2DIndex(i + paddingSize, j + paddingSize, newShape),))
        val.append(sf[i, j])

field = interpolateScalarField(newShape, P, verts,
                               weights, val,
                               [], [], [],
                               fixBoundary=True, boundaryValue=paddingBoundaryVal,
                               dType=np.float64, )

padded = field.reshape(*newShape)

assert np.max(np.abs(padded[paddingSize:-paddingSize, paddingSize:-paddingSize] - smoothed)) < 1e-4

# np.save(outfile + ".npy", smoothed)
# imageToVTK(outfile, pointData={"speed": speed_pts, "speed_smoothed": smoothed_pts})
