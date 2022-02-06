import pyvista as pv
import numpy as np
import itertools, glob, os, tqdm
from pathlib import Path
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

def drawMergeTree(pts, edges, outFile):

    # Draw the tree
    polyData = vtk.vtkPolyData()
    ptsVtk = vtk.vtkPoints()

    for i in range(len(pts)):
        ptsVtk.InsertNextPoint(pts[i])
    polyData.SetPoints(ptsVtk)
    lines = vtk.vtkCellArray()

    for i in range(len(edges)):
        line = vtk.vtkLine()

        line.GetPointIds().SetId(0, edges[i][0])  # the second 0 is the index of the Origin in the vtkPoints
        line.GetPointIds().SetId(1, edges[i][1])  # the second 1 is the index of P0 in the vtkPoints
        lines.InsertNextCell(line)

    polyData.SetLines(lines)
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polyData)
    writer.SetFileName(join(outFile))
    writer.Update()


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

def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)

def contourLineConstraint(O, r, disThreshold=0.5):
    coords = []
    for i,j in itertools.product(range(200), range(200)):
        if i <= O[0] + r and i >= O[0] - r and j <= O[1] + r and j >= O[1] - r:
            coords .append(np.array([i,j]))
    coords = np.array(coords)
    # determine if pixel's distance is less than 0.5
    ## x direction
    toCircleDisXDir = np.abs(np.sqrt(r**2 - (O[1] - coords[:, 1])**2) - np.abs(O[0] - coords[:, 0]))
    toCircleDisYDir = np.abs(np.sqrt(r**2 - (O[0] - coords[:, 0])**2) - np.abs(O[1] - coords[:, 1]))

    pixelsOnCircle = np.union1d(np.where(toCircleDisXDir<disThreshold)[0], np.where(toCircleDisYDir<disThreshold)[0])
    coordsOnCirlcle = coords[pixelsOnCircle]

    return coordsOnCirlcle

def contourLineConstraintOnInputCoords(coords, O, r, disThreshold=0.5):
    # determine if pixel's distance is less than 0.5
    ## x direction
    toCircleDisXDir = np.abs(np.sqrt(r**2 - (O[1] - coords[:, 1])**2) - np.abs(O[0] - coords[:, 0]))
    toCircleDisYDir = np.abs(np.sqrt(r**2 - (O[0] - coords[:, 0])**2) - np.abs(O[1] - coords[:, 1]))

    pixelsOnCircle = np.union1d(np.where(toCircleDisXDir<disThreshold)[0], np.where(toCircleDisYDir<disThreshold)[0])

    return pixelsOnCircle

def contourLineInterpolationConstraints(O, r, gridSize = (200, 200)):

    constraintIds = []
    constraintWeights = []

    # intersection with axis 0
    for i in range(int(np.ceil(O[0]-r)), int(np.floor(O[0]+r))):
        dy = np.sqrt(r**2 - (i-O[0])**2)
        if dy > 1e-5:
            j = O[1] + dy
            constraintIds.append([flatten2DIndex(i, np.floor(j), gridSize=gridSize), flatten2DIndex(i, np.ceil(j), gridSize=gridSize)])
            constraintWeights.append([1-(j-np.floor(j)), 1-(np.ceil(j)-j)])

            j = O[1] - dy
            constraintIds.append([flatten2DIndex(i, np.floor(j), gridSize=gridSize), flatten2DIndex(i, np.ceil(j), gridSize=gridSize)])
            constraintWeights.append([1-(j-np.floor(j)), 1-(np.ceil(j)-j)])
        else:
            j = O[1]
            constraintIds.append(
                [flatten2DIndex(i, np.floor(j), gridSize=gridSize), flatten2DIndex(i, np.ceil(j), gridSize=gridSize)])
            constraintWeights.append([1 - (j - np.floor(j)), 1 - (np.ceil(j)-j)])

    # intersection with axis 0
    for j in  range(int(np.ceil(O[1]-r)), int(np.floor(O[1]+r))):
        dy = np.sqrt(r ** 2 - (j - O[1]) ** 2)
        if dy > 1e-5:
            i = O[0] + dy
            constraintIds.append([flatten2DIndex(np.floor(i), j, gridSize=gridSize),
                                  flatten2DIndex(np.ceil(i), j, gridSize=gridSize)])
            constraintWeights.append([1 - (i - np.floor(i)), 1 - (np.ceil(i)-i)])

            i = O[0] - dy
            constraintIds.append([flatten2DIndex(np.floor(i), j, gridSize=gridSize),
                                  flatten2DIndex(np.ceil(i), j, gridSize=gridSize)])
            constraintWeights.append([1 - (i - np.floor(i)), 1 - (np.ceil(i)-i)])
        else:
            i = O[0]
            constraintIds.append([flatten2DIndex(np.floor(i), j, gridSize=gridSize),
                                  flatten2DIndex(np.ceil(i), j, gridSize=gridSize)])
            constraintWeights.append([1 - (i - np.floor(i)), 1 - (np.ceil(i) - i)])

    return constraintIds, constraintWeights

def obj2vtkFolder(inObjFolder, inFileExt='obj', outVtkFolder=None, processInterval=[], addFaces=False,
                  addABeforeName=True, faceMesh=None, ):
    # addFaces = True
    if outVtkFolder is None:
        outVtkFolder = inObjFolder + r'\vtk'

    objFiles = glob.glob(inObjFolder + r'\*.' + inFileExt)


    meshWithFaces = None

    os.makedirs(outVtkFolder, exist_ok=True)
    if len(processInterval) == 2:
        objFiles = objFiles[processInterval[0]: processInterval[1]]
    for f in tqdm.tqdm(objFiles, desc=inFileExt + " to vtk"):
        fp = Path(f)

        mesh = pv.read(f)
        if addFaces and meshWithFaces is not None:
            mesh.faces = meshWithFaces.faces
        # else:
        #     mesh.faces = np.empty((0,), dtype=np.int32)
        if addABeforeName:
            outName = outVtkFolder + r'\\A' + fp.stem + '.vtk'
        else:
            outName = outVtkFolder + r'\\' + fp.stem + '.vtk'

        mesh.point_arrays['Height'] = mesh.points[:, 2]

        mesh.save(outName)


def genLaplacianMat(gridSize):
    numViables = gridSize[0] * gridSize[1]
    A = np.zeros((numViables, numViables,), dtype=np.float32)
    # A = np.zeros((numViables, numViables, ), dtype=np.float64)

    for i, j in itertools.product(range(1, gridSize[0]-1), range(1, gridSize[1]-1)):
        # y direction second derivatives
        id_ij = flatten2DIndex(i, j, gridSize)
        id_im1j = flatten2DIndex(i - 1, j, gridSize)
        id_ip1j = flatten2DIndex(i + 1, j, gridSize)

        A[id_ij, id_ij] += 4
        A[id_im1j, id_im1j] += 1
        A[id_ip1j, id_ip1j] += 1

        A[id_ij, id_ip1j] += -2
        A[id_ip1j, id_ij] += -2

        A[id_ij, id_im1j] += -2
        A[id_im1j, id_ij] += -2

        A[id_ip1j, id_im1j] += 1
        A[id_im1j, id_ip1j] += 1

        # x direction second derivatives
        id_ijm1 = flatten2DIndex(i, j - 1, gridSize)
        id_ijp1 = flatten2DIndex(i, j + 1, gridSize)

        A[id_ij, id_ij] += 4
        A[id_ijm1, id_ijm1] += 1
        A[id_ijp1, id_ijp1] += 1

        A[id_ij, id_ijp1] += -2
        A[id_ijp1, id_ij] += -2

        A[id_ij, id_ijm1] += -2
        A[id_ijm1, id_ij] += -2

        A[id_ijp1, id_ijm1] += 1
        A[id_ijm1, id_ijp1] += 1

    return A


def check2DCoordValidility(i, j, gridSize):

    return i >= 0 and i < gridSize[0] and j >= 0 and j < gridSize[1]

def findContourLineHeight(nodes, fieldSeg, gridSize, directions):
    contourLineHeight = {} # key (a, b): a: higher segment, c lower segment
    iMeshVerts = matchTreeToGrid(nodes.points, fieldSeg.points)

    for iNode in range(nodes.points.shape[0]):

        if nodes['CriticalType'][iNode] != 2:
            continue
        # how to find the height of a contour line?
        # the height of a contour line is the height of a saddle point where to contour line meets
        # at the neighborhood of each saddle point, there will be three types of nodes: two higher segment and one lower segment
        # we will need to find the two higher seg, say (a,b) and the lower seg c
        # the height of contour line a-c, b-c will be the height (scalar) of this saddle point

        # match node to grid mesh
        iMeshVert = iMeshVerts[iNode]
        i, j = to2DIndex(iMeshVert, gridSize)

        segsInNeighborhood = []
        heights = []

        for d in directions:
            if check2DCoordValidility(i+d[0], j+d[1], gridSize, ):
                neiVertId = flatten2DIndex(i+d[0], j+d[1], gridSize)
                if fieldSeg['SegmentationId'][neiVertId] not in segsInNeighborhood:
                    segsInNeighborhood.append(fieldSeg['SegmentationId'][neiVertId])
                    heights.append(fieldSeg['Height'][neiVertId])

        assert len(segsInNeighborhood) == 3
        sortedId = np.argsort(heights)
        contourLineHeight[(segsInNeighborhood[sortedId[2]], segsInNeighborhood[sortedId[0]])] = nodes['Scalar'][iNode]
        contourLineHeight[(segsInNeighborhood[sortedId[1]], segsInNeighborhood[sortedId[0]])] = nodes['Scalar'][iNode]

    return contourLineHeight