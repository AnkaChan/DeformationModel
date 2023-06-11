import pyvista as pv
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from pyevtk.hl import imageToVTK

input_dir = "./VortexSlice/"
output_dir = "./VortexSlice/monoMesh_8_10_smoothed_padded/"
inputStartFile = "monoMesh_010_smoothed_padded.npy"
outfile = input_dir + "monoMesh_010_flipped"
# gridSize = (192, 64)

fig = plt.figure()
# input = pv.read(input_dir + inputStartFile)
input_image = np.load(input_dir + inputStartFile)
gridSize = input_image.shape

# input.plot()
# input_image = np.array(input['Scalars_']).reshape(gridSize[1], gridSize[0]).T
input_image_flipped = -1*input_image
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.imshow(input_image)
ax2.imshow(input_image_flipped)
plt.show()

input_pts = input_image.reshape(gridSize[0], gridSize[1], 1)
flipped_pts = input_image_flipped.reshape(gridSize[0], gridSize[1], 1)
imageToVTK(outfile, pointData={"speed": input_pts, "speed_neg": flipped_pts})