import pyvista as pv
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from pyevtk.hl import imageToVTK

input_dir = "./VortexSlice/"
# output_dir = "./VortexSlice/monoMesh_8_10/"
inputStartFile = "monoMesh_008.vti"
outfile = input_dir + "monoMesh_008_smoothed"
gridSize = (192, 64)

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
speed_pts = input_image.reshape(gridSize[0], gridSize[1], 1)
smoothed_pts = smoothed.reshape(gridSize[0], gridSize[1], 1)
print(outfile)
np.save(outfile + ".npy", smoothed)
# imageToVTK(outfile, pointData={"speed": speed_pts, "speed_smoothed": smoothed_pts})
