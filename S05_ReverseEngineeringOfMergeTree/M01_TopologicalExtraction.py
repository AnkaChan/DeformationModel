import pyvista as pv
import numpy as np
import itertools, glob, os, tqdm
from pathlib import Path
import cv2
from os.path import join
from scipy import sparse
from vtk.util.numpy_support import vtk_to_numpy
import vtk

from qpsolvers import solve_qp, available_solvers
from scipy import sparse
print('Avaliable qp solvers: ', available_solvers)
def writeOBj(outObj, N, X, Y, Z):
    '''
    Write a scalar field to a grid mesh
    '''
    file  = open(outObj, 'w')
    for i, j in itertools.product(range(N), range(N)):
        file.write('v %f %f %f\n' %( X[i, j],  Y[i,j], Z[i,j] ))

    for i, j in itertools.product(range(0, N-1), range(1, N)):
        vId = j + i *N
        file.write('f %d %d %d\n' %(vId, vId+1,  vId+N+1, ))
        file.write('f %d %d %d\n' %(vId, vId+N+1,  vId+N, ))

def flatten2DIndex(i,j, gridSize, major='row'):
    '''
    i: row
    j: col
    gridSize: (rows, cols)
    '''
    if major == 'col':
        return gridSize[0] * j + i
    elif major == 'row':
        return gridSize[1] * i + j

def to2DIndex(id, gridSize, major='row'):
    '''
    i: row
    j: col
    '''
    if major == 'col':
        return (int(id % gridSize[0]), int((np.floor(id / gridSize[0]))))
    elif major == 'row':
        return  (int(np.floor(id / gridSize[1])), int(id % gridSize[1]))
def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

def matchTreeToGrid(treeNodes, flattenGrid):
    numNodes = treeNodes.shape[0]

    matchedVIds = []
    matchDis = []
    for iNode in range(numNodes):
        diff = flattenGrid - treeNodes[iNode, :]
        dis = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)

        matchedVIds.append(np.argmin(dis))
        matchDis.append(dis[matchedVIds[-1]])
        print("Min distance: ", dis[matchedVIds[-1]])

    return matchedVIds

def getVTIData(vtiImage, arrayName, reshape2D=True):
    '''
    vtiImage: vtkImageData
    reshape2D: whether reshape to a 2D array, otherwise flatten
    return: flatten array with row major if not reshape2D
    '''
    dims = vtiImage.GetDimensions()
    # dims = (width, height) = (cols, rows)
    dimsRC = (dims[1], dims[0])

    scalar = vtk_to_numpy(vtiImage.GetPointData().GetArray(arrayName))

    if reshape2D:
        img = scalar.reshape(dimsRC)
        # ‘C’ means to read / write the elements using C-like index order,
        # with the last axis index changing fastest, back to the first axis index changing slowest.
        # such scalar is row major
        return img, dimsRC
    else:
        return scalar, dimsRC

def readVTIData(vtiImageFile, arrayName, reshape2D=True):
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(vtiImageFile)

    reader.Update()
    image = reader.GetOutput()

    return getVTIData(image, arrayName, reshape2D=reshape2D)

def findContourLineConstraints(fieldSeg, gridSize, contourLineHeight, edges=None):
    '''
    fieldSeg: flatten segmentation Ids, row major
    gridSize: (rows, cols)
    contourLineHeight: a dict contains contour line height between two segments. (seg1, seg2): height. Note that seg1 >= seg2
    '''
    # don't consider boundary edges
    directionsForEdges = [
        (1, 0),
        (0, 1),
    ]
    if edges is None:
        edges = []
        for i,j in tqdm.tqdm(itertools.product(range(gridSize[0]), range(gridSize[1]))):
            for d in directionsForEdges:
                if i + d[0] >= 0 and i + d[0] < gridSize[0] and j + d[1] >= 0 and j + d[1] < gridSize[1]:
                    edges.append([flatten2DIndex(i, j, gridSize=gridSize), flatten2DIndex(i+d[0], j+d[1], gridSize=gridSize)])

    contourEdges = []
    contourLineConstraintWeight = []
    contourLineConstraintHeight = []
    for i, edge in tqdm.tqdm(enumerate(edges)):
        if fieldSeg["SegmentationId"][edge[0]] != fieldSeg["SegmentationId"][edge[1]]:

            contourEdges.append(edge)
            h1 = fieldSeg["Height"][edge[0]]
            h2 = fieldSeg["Height"][edge[1]]
            if h1 > h2:
                if contourLineHeight.get(
                        (fieldSeg["SegmentationId"][edge[0]], fieldSeg["SegmentationId"][edge[1]])) is None:
                    continue
                contourHeight = contourLineHeight[fieldSeg["SegmentationId"][edge[0]], fieldSeg["SegmentationId"][edge[1]]]
            else:
                if contourLineHeight.get(
                        (fieldSeg["SegmentationId"][edge[1]], fieldSeg["SegmentationId"][edge[0]])) is None:
                    continue
                contourHeight = contourLineHeight[
                    fieldSeg["SegmentationId"][edge[1]], fieldSeg["SegmentationId"][edge[0]]]

            w1 = (contourHeight - h2) / (h1 - h2)
            assert 0 <= w1 <=1
            w2 = 1 - w1
            contourLineConstraintWeight.append([w1, w2])
            contourLineConstraintHeight.append(contourHeight)

    return contourEdges, contourLineConstraintWeight, contourLineConstraintHeight