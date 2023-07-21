import pyvista as pv

inFile = r'contour_14.vtp'

data = pv.PolyData(inFile)

print(data)