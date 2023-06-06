import array

import pyvista as pv
import numpy as np
import itertools, glob, os, tqdm
from pathlib import Path
from os.path import join
from scipy import sparse
import matplotlib.pyplot as plt

from vtk.util.numpy_support import vtk_to_numpy
import vtk

from qpsolvers import solve_qp, available_solvers
from scipy import sparse
from matplotlib.collections import LineCollection

from shapely.geometry.polygon import Polygon, Point
print('Avaliable qp solvers: ', available_solvers)
# def writeOBj(outObj, N, X, Y, Z):
#     '''
#     Write a scalar field to a grid mesh
#     '''
#     file  = open(outObj, 'w')
#     for i, j in itertools.product(range(N), range(N)):
#         file.write('v %f %f %f\n' %( X[i, j],  Y[i,j], Z[i,j] ))
#
#     for i, j in itertools.product(range(0, N-1), range(1, N)):
#         vId = j + i *N
#         file.write('f %d %d %d\n' %(vId, vId+1,  vId+N+1, ))
#         file.write('f %d %d %d\n' %(vId, vId+N+1,  vId+N, ))


class TreeNode:
    def __init__(s,):
        s.id = None
        s.criticalType = None
        s.position = [] # for visulization
        s.scalar = None
        s.posInField = None # position in terms of row, col at the scalar field

        s.upNodes = []
        s.downNodes = []
        s.upEdges = []
        s.downEdges = []

        # for intermediate tree
        s.tree0Corr = -1
        s.tree1Corr = -1

        #
        s.type = None
        # saddle point only, value can be partial or complete
        # partial: one contourline is preserved
        # complete: no contourlien is preserved
        s.emergeVanishType = None
        s.preservingContourId = None

class Edge:
    def __init__(s,):
        s.id = None
        s.segmentId = None
        s.upNode = None
        s.downNode = None
        s.nodes = []

        # for intermediate tree
        s.segmentIdTree0 = -1
        s.segmentIdTree1 = -1


class ContourLine:
    def __init__(s, gridSize):
        s.saddleAllContourEdges = None
        s.saddleAllContourWeights = None
        s.saddleAllContourHeights = None
        s.embracingHigherNodeId = None
        s.contourLineParameters = []
        s.shapelyGeometry = None
        s.relativeTranslations = None

        s.gridSize = gridSize

        s.allNodes = None

        s.upSegment = None
        s.downSegment = None

        s.startState = None
        s.endState = None

        s.type = None

        s.scalar = None

    def computeRelativeTranslation(s):
        # contour line always start and end with saddle
        s.relativeTranslations = []

        saddlePoint = s.allNodes[0]
        for iNdoe in range(len(s.allNodes)):
            diff = s.allNodes[iNdoe] - saddlePoint
            s.relativeTranslations.append(diff)
        s.relativeTranslations = np.array(s.relativeTranslations)

    def getRelativeTranslation(s, t):
        ts = s.contourLineParameters

        if t >= 1:
            t = t-1

        diffa = ts - t
        diffa[diffa >0] = -2
        a = np.argmax(diffa)

        diffb = ts - t
        diffb[diffb <0] = 2
        b = np.argmin(diffb)
        # a = np.max(np.where(t >= ts))
        # b = np.min(np.where(t <= ts))

        if b is None:
            b = a
        if a == b:
            return s.relativeTranslations[a]
        # assert a == b - 1
        controlPoints = s.relativeTranslations
        w = (t - ts[a]) / (ts[b] - ts[a])

        return controlPoints[a] * (1-w) + controlPoints[b] * w

    def getSaddle(s):
        return s.allNodes[0]

    def numVertices(s):
        return len(s.saddleAllContourEdges)

    def initializeShapely(s):
        s.shapelyGeometry = Polygon(s.allNodes)

    def computeNode(s, iNode, ):
        p1 = np.array(to2DIndex(s.saddleAllContourEdges[iNode][0], s.gridSize))
        p2 = np.array(to2DIndex(s.saddleAllContourEdges[iNode][1], s.gridSize))
        p = p1 * s.saddleAllContourWeights[iNode][0] \
            + p2 * s.saddleAllContourWeights[iNode][1]
        return p

    def getNode(s, iNode):
        if s.allNodes is None:
            s.initializeAllNodes()
        return s.allNodes[iNode]

    def getAllNodes(s):
        if s.allNodes is None:
            s.initializeAllNodes()
        return s.allNodes

    def initializeAllNodes(s):
        allPts = []
        n = s.numVertices()

        for iNode in range(n):
            allPts.append(s.computeNode(iNode))

        s.allNodes = np.array(allPts)

    def intialize(s):
        s.initializeAllNodes()
        s.initializeShapely()

    def getEdge(s,  iEdge):
        if iEdge == len(s.saddleAllContourEdges) - 1:
            nextNode = 0
        else:
            nextNode = iEdge + 1

        return  s.getNode(nextNode) - s.getNode(iEdge)

    def checkCCW(s):
        totalOrientedArea = 0
        for i in range(1, len(s.saddleAllContourEdges)-1):
            startingP = s.getNode(0)
            edgeP1 = s.getNode(i)
            nextEdge = s.getEdge(i)

            p0p1 = startingP - edgeP1

            totalOrientedArea = totalOrientedArea + nextEdge[0] * p0p1[1] -  nextEdge[1] * p0p1[0]

        if totalOrientedArea > 0:
            return True
        else:
            return False

    def parameterize(s, parameterStartNodeId=0, repermute=False):
        totalLength = 0
        s.contourLineParameters = np.zeros((s.numVertices(),))
        for iVert in range(s.numVertices()-1):
            vertId = iVert + parameterStartNodeId
            if vertId >= s.numVertices():
                vertId = vertId - s.numVertices()

            edgeLength = np.linalg.norm(s.getEdge(vertId))
            s.contourLineParameters[vertId] = (totalLength + edgeLength)
            totalLength += edgeLength

        totalLength += np.linalg.norm(s.getEdge(s.numVertices()-1))
        s.contourLineParameters = s.contourLineParameters/totalLength

        if parameterStartNodeId != 0 and repermute:
            # sort the nodes to make sure their parameters are ascending
            sortedId = np.argsort(s.contourLineParameters)
            s.contourLineParameters = s.contourLineParameters[sortedId]
            s.saddleAllContourEdges = [s.saddleAllContourEdges[id] for id in sortedId]
            s.saddleAllContourWeights = [s.saddleAllContourWeights[id] for id in sortedId]
            s.saddleAllContourHeights =[s.saddleAllContourHeights[id] for id in sortedId]

            # s.saddleAllContourEdges = s.saddleAllContourEdges[sortedId]
            # s.saddleAllContourWeights = s.saddleAllContourWeights[sortedId]
            # s.saddleAllContourHeights = s.saddleAllContourHeights[sortedId]
            # s.embracingHigherNodeId = s.embracingHigherNodeId[sortedId]
            s.initializeAllNodes()

    def getPosition(s, t):
        ts = s.contourLineParameters

        if t >= 1:
            t = t-1

        diffa = ts - t
        diffa[diffa >0] = -2
        a = np.argmax(diffa)

        diffb = ts - t
        diffb[diffb <0] = 2
        b = np.argmin(diffb)
        # a = np.max(np.where(t >= ts))
        # b = np.min(np.where(t <= ts))

        if b is None:
            b = a
        if a == b:
            return s.getNode(a)
        # assert a == b - 1
        controlPoints = s.allNodes
        w = (t - ts[a]) / (ts[b] - ts[a])

        return controlPoints[a] * (1-w) + controlPoints[b] * w

    def reverseOrientation(s):
        s.saddleAllContourEdges.reverse()
        s.saddleAllContourWeights.reverse()
        s.saddleAllContourHeights.reverse()
        s.contourLineParameters.reverse()


    def includePoint(s, p):
        n = s.numVertices()
        inside = False

        p1 = s.getNode(0)
        p1x, p1y = p1[0], p1[1]
        for i in range(n + 1):
            p2x, p2y =s.getNode( i % n)
            if p[1] > min(p1y, p2y):
                if p[1] <= max(p1y, p2y):
                    if p[0] <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (p[1] - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or p[0] <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y

        return inside

class ContourConstraint:
    def __init__(s, grideSize):
        s.contourLines = []
        s.gridSize = grideSize

    def getContour(s, iContour):
        return s.contourLines[iContour]


    def getNode(s, iContour, iNode):
        return s.getContour(iContour).getNode(iNode)

    def getAllNodes(s,iContour):
        return s.getContour(iContour).getAllNodes()

    def numContours(s):
        return len(s.contourLines)

    def numVertices(s, iContour):
        return len(s.getContour(iContour).saddleAllContourEdges[iContour])

    def computeNode(s, iContour, iNode, ):
        p1 = np.array(to2DIndex(s.getContour(iContour).saddleAllContourEdges[iContour][iNode][0], s.gridSize))
        p2 = np.array(to2DIndex(s.getContour(iContour).saddleAllContourEdges[iContour][iNode][1], s.gridSize))
        p = p1 * s.getContour(iContour).saddleAllContourWeights[iContour][iNode][0] \
            + p2 * s.getContour(iContour).saddleAllContourWeights[iContour][iNode][1]
        return p

    def getEdge(s, iContour,  iEdge):
        if iEdge == len(s.getContour(iContour).saddleAllContourEdges) - 1:
            nextNode = 0
        else:
            nextNode = iEdge + 1

        return  s.getNode(iContour, nextNode) - s.getNode(iContour, iEdge)



    def includePoint(s, iContour, p):
        n = s.numVertices(iContour)
        inside = False

        p1 = s.getNode(iContour, 0)
        p1x, p1y = p1[0], p1[1]
        for i in range(n + 1):
            p2x, p2y =s.getNode(iContour, i % n)
            if p[1] > min(p1y, p2y):
                if p[1] <= max(p1y, p2y):
                    if p[0] <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (p[1] - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or p[0] <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y

        return inside


    def parameterize(s):
        for iContour in range(s.numContours()):
            s.getContour(iContour).contourLineParameters.append([])
            totalLength = 0

            for iVert in range(s.numVertices(iContour)):
                s.getContour(iContour).contourLineParameters[-1].append(np.linalg.norm(s.getNode(iContour, iVert)))
                totalLength += s.getContour(iContour).contourLineParameters[-1][-1]

            s.getContour(iContour).contourLineParameters[-1] = np.array(s.getContour(iContour).contourLineParameters[-1]/totalLength)

    def getPosition(s, iContour, t):
        ts = s.getContour(iContour).contourLineParameters[iContour]

        a = np.max(np.where(t >= ts))
        b = np.min(np.where(t <= ts))

        if b is None:
            b = a + 1
        if a == b:
            return s.getNode(iContour, a)
        assert a == b - 1
        controlPoints = s.getContour(iContour).allNodes[iContour]

        w = (t - ts[a]) / (ts[b] - ts[a])

        return controlPoints[a] * (1-w) + controlPoints[b] * w

    def addContour(s, newContour):
        s.contourLines.append(newContour)

class Tree:
    def __init__(s, nodesData=None, edgesData=None, segmentationData=None, gridSize=None, splitTree=False, segmentationDataScalarName="Height"):
        s.nodes = []
        s.edges = []

        s.isSplitTree = splitTree

        s.rootNodeId = None

        if splitTree:
            actualUp = "down"
            actualDown = "up"
            if nodesData is not None:
                nodesData['Scalar'] = - nodesData['Scalar']
            if segmentationData is not None:
                segmentationData[segmentationDataScalarName] = -segmentationData[segmentationDataScalarName]

                # scalars = np.array(segmentationData[segmentationDataScalarName])
                # plt.imshow(scalars.reshape(gridSize))
                # plt.show()
        else:
            actualUp = "up"
            actualDown = "down"


        s.nodeToField = []
        s.contourLineHeight = {}  # key (a, b): a: higher segment, b lower segment

        s.nodesData = nodesData
        s.edgesData = edgesData
        s.segmentationData = segmentationData
        s.segmentationDataScalarName = segmentationDataScalarName
        s.saddleContours = {}

        s.saddleTypeId = 2

        # convert split tree to merge tree
        # if splitTree:
        #     s.nodesData['Scalar'] = -s.nodesData['Scalar']
        #     s.segmentationData[s.segmentationDataScalarName] = -s.segmentationData[s.segmentationDataScalarName]

        if gridSize is None:
            s.gridSize = ()
        else:
            s.gridSize = gridSize


        if nodesData is not None and edgesData is not None:

            for iNode in range(nodesData.points.shape[0]):
                # non-saddle node
                newNode = TreeNode()
                newNode.id = iNode
                newNode.criticalType = nodesData['CriticalType'][iNode]
                newNode.scalar = nodesData['Scalar'][iNode]
                newNode.position = nodesData.points[iNode]
                s.nodes.append(newNode)

            for iEdge in range(edgesData["upNodeId"].shape[0]):
                newEdge = Edge()
                newEdge.id = iEdge
                newEdge.segmentId = edgesData["SegmentationId"][iEdge]
                newEdge.nodes = [edgesData[actualUp + "NodeId"][iEdge], edgesData[actualDown+"NodeId"][iEdge],]
                newEdge.upNode = edgesData[actualUp + "NodeId"][iEdge]
                newEdge.downNode = edgesData[actualDown+"NodeId"][iEdge]
                s.edges.append(newEdge)
                # initialize the connectivity infos for up node

                upNode = s.nodes[newEdge.upNode]
                upNode.downEdges.append(iEdge)
                upNode.downNodes.append(edgesData[actualDown+"NodeId"][iEdge])

                downNode = s.nodes[newEdge.downNode]
                downNode.upNodes.append(edgesData[actualUp+"NodeId"][iEdge])
                downNode.upEdges.append(iEdge)

                # if edgesData["downNodeId"][iEdge] == iNode:
                #     upSegs.append(edgesData["SegmentationId"][iEdge])
                # if edgesData["upNodeId"][iEdge] == iNode:
                #     downsegs.append(edgesData["SegmentationId"][iEdge])

        if segmentationData is not None:
            s.matchTreeToGrid()

    def numNodes(s):
        return len(s.nodes)

    def node(s, iNode):
        return s.nodes[iNode]

    def numEdges(s):
        return len(s.edges)


    def initFrom(s, nodes, edges):
        s.edges = edges
        s.nodes = nodes

    def load(s, nodeFile, edgeFile, segmentationFile, gridSize=None, splitTree=False, segmentationDataScalarName="Height"):
        fieldSeg = pv.read(segmentationFile)
        nodes = pv.read(nodeFile)
        edges = pv.read(edgeFile)
        s.__init__(nodes, edges, fieldSeg, gridSize=gridSize, splitTree=splitTree, segmentationDataScalarName=segmentationDataScalarName)

    def matchTreeToGrid(s ):
        numNodes = len(s.nodes)

        matchDis = []
        for iNode in range(numNodes):
            diff = s.segmentationData.points - s.nodes[iNode].position
            dis = np.sqrt(diff[:, 0] ** 2 + diff[:, 1] ** 2 + diff[:, 2] ** 2)

            s.nodeToField.append(np.argmin(dis))
            matchDis.append(dis[s.nodeToField[-1]])
            # print("Min distance: ", dis[matchedVIds[-1]])
            assert dis[s.nodeToField[-1]] < 1e-6

            s.nodes[iNode].posInField = np.array(to2DIndex(s.nodeToField[-1], s.gridSize))


    def findContourLineHeight(s, ):
        """
        In the merge tree, nodes corresponds to critical points and the edges corresponds to the segments
        """

        for iNode in range(s.nodesData.points.shape[0]):
            # non-saddle node
            if s.nodesData['CriticalType'][iNode] != s.saddleTypeId:
                continue
            # saddle node
            # iterate all edges to find two nodes above the node
            upSegs = []
            downsegs = []

            for iEdge in range(len(s.edges)):
                if s.edges[iEdge].downNode == iNode:
                    upSegs.append(s.edges[iEdge].segmentId)
                if s.edges[iEdge].upNode == iNode:
                    downsegs.append(s.edges[iEdge].segmentId)

            assert len(upSegs) == 2
            assert len(downsegs) == 1

            s.contourLineHeight[(upSegs[0], downsegs[0])] = s.nodesData["Scalar"][iNode]
            s.contourLineHeight[(upSegs[1], downsegs[0])] = s.nodesData["Scalar"][iNode]


    def extractContourLineConstraints(s, edges=None):
        '''
        fieldSeg: flatten segmentation Ids, row major
        gridSize: (rows, cols)
        contourLineHeight: a dict contains contour line height between two segments. (seg1, seg2): height. Note that seg1 >= seg2
        edges: candidate edges to be checked, if it's set none, it will check all the edges

        '''

        s.findContourLineHeight()
        # don't consider boundary edges
        directionsForEdges = [
            (1, 0),
            (0, 1),
        ]
        if edges is None:
            edges = []
            for i,j in tqdm.tqdm(itertools.product(range(s.gridSize[0]), range(s.gridSize[1])), desc="Generating candidate edges."):
                for d in directionsForEdges:
                    if i + d[0] >= 0 and i + d[0] < s.gridSize[0] and j + d[1] >= 0 and j + d[1] < s.gridSize[1]:
                        assert flatten2DIndex(i, j, gridSize=s.gridSize) == flatten2DIndex(*to2DIndex(flatten2DIndex(i, j, gridSize=s.gridSize), s.gridSize), s.gridSize)
                        edges.append((flatten2DIndex(i, j, gridSize=s.gridSize), flatten2DIndex(i+d[0], j+d[1], gridSize=s.gridSize)))

        s.contourIntersectingEdges = []
        s.contourLineConstraintWeight = []
        s.contourLineConstraintHeight = []
        s.contourLineConstraintUpDownSegments = []
        for i, edge in tqdm.tqdm(enumerate(edges), desc="Examining all edges for counter line intersection."):
            if s.segmentationData["SegmentationId"][edge[0]] != s.segmentationData["SegmentationId"][edge[1]]:
                s.contourIntersectingEdges.append(edge)
                h1 = s.segmentationData[s.segmentationDataScalarName][edge[0]]
                h2 = s.segmentationData[s.segmentationDataScalarName][edge[1]]
                if h1 > h2:
                    contourNeiSegs = (s.segmentationData["SegmentationId"][edge[0]], s.segmentationData["SegmentationId"][edge[1]])
                    if s.contourLineHeight.get(contourNeiSegs) is None:
                        s.contourIntersectingEdges.pop()
                        print("Warining! No contour line height defined between: ", contourNeiSegs)
                        continue
                    contourHeight = s.contourLineHeight[s.segmentationData["SegmentationId"][edge[0]], s.segmentationData["SegmentationId"][edge[1]]]
                    contoutUpDownSegments = (s.segmentationData["SegmentationId"][edge[0]], s.segmentationData["SegmentationId"][edge[1]])
                else:
                    contourNeiSegs = (s.segmentationData["SegmentationId"][edge[1]], s.segmentationData["SegmentationId"][edge[0]])
                    if s.contourLineHeight.get(contourNeiSegs) is None:
                        s.contourIntersectingEdges.pop()
                        print("Warining! No contour line height defined between: ", contourNeiSegs)

                        continue
                    contourHeight = s.contourLineHeight[
                        s.segmentationData["SegmentationId"][edge[1]], s.segmentationData["SegmentationId"][edge[0]]]
                    contoutUpDownSegments = (s.segmentationData["SegmentationId"][edge[1]], s.segmentationData["SegmentationId"][edge[0]])

                # check if h1 == h2
                if h1 == h2:
                    w1 = 0.5
                else:
                    w1 = (contourHeight - h2) / (h1 - h2)

                assert 0 <= w1 <=1
                w2 = 1 - w1

                s.contourLineConstraintWeight.append((w1, w2))
                s.contourLineConstraintHeight.append(contourHeight)
                s.contourLineConstraintUpDownSegments.append(contoutUpDownSegments)

            if len(s.contourIntersectingEdges) != len(s.contourLineConstraintWeight):
                break

    def extractContourLineConstraintsNew(s, edges=None):
        '''
        fieldSeg: flatten segmentation Ids, row major
        gridSize: (rows, cols)
        contourLineHeight: a dict contains contour line height between two segments. (seg1, seg2): height. Note that seg1 >= seg2
        edges: candidate edges to be checked, if it's set none, it will check all the edges

        '''

        s.findContourLineHeight()
        # don't consider boundary edges
        directionsForEdges = [
            (1, 0),
            (0, 1),
        ]
        # if edges is None:
        #     edges = []
        #     for i, j in tqdm.tqdm(itertools.product(range(s.gridSize[0]), range(s.gridSize[1])),
        #                           desc="Generating candidate edges."):
        #         for d in directionsForEdges:
        #             if i + d[0] >= 0 and i + d[0] < s.gridSize[0] and j + d[1] >= 0 and j + d[1] < s.gridSize[1]:
        #                 assert flatten2DIndex(i, j, gridSize=s.gridSize) == flatten2DIndex(
        #                     *to2DIndex(flatten2DIndex(i, j, gridSize=s.gridSize), s.gridSize), s.gridSize)
        #                 edges.append((flatten2DIndex(i, j, gridSize=s.gridSize),
        #                               flatten2DIndex(i + d[0], j + d[1], gridSize=s.gridSize)))
        segmentedField2D =s.segmentationData["SegmentationId"].reshape(s.gridSize)

        allContourEdges = []
        for d in directionsForEdges:
            s1 = segmentedField2D[:segmentedField2D.shape[0]-d[0], :segmentedField2D.shape[1]-d[1]]
            s2 = segmentedField2D[d[0]:, d[1]:]

            contouredges2D = np.where(s1 != s2)
            contourEdgesV1Flatten = flatten2DIndex(contouredges2D[0], contouredges2D[1], s.gridSize)
            contourEdgesV2Flatten = flatten2DIndex(contouredges2D[0] + d[0], contouredges2D[1] + d[1], s.gridSize)

            allContourEdges.extend([(contourEdgesV1Flatten[iE], contourEdgesV2Flatten[iE]) for iE in range(contourEdgesV1Flatten.shape[0])])

        s.contourIntersectingEdges = []
        s.contourLineConstraintWeight = []
        s.contourLineConstraintHeight = []
        s.contourLineConstraintUpDownSegments = []

        for iContourEdge in range(len(allContourEdges)):
            edge = allContourEdges[iContourEdge]
            if s.segmentationData["SegmentationId"][edge[0]] != s.segmentationData["SegmentationId"][edge[1]]:
                s.contourIntersectingEdges.append(edge)
                h1 = s.segmentationData[s.segmentationDataScalarName][edge[0]]
                h2 = s.segmentationData[s.segmentationDataScalarName][edge[1]]
                if h1 > h2:
                    contourNeiSegs = (
                    s.segmentationData["SegmentationId"][edge[0]], s.segmentationData["SegmentationId"][edge[1]])
                    if s.contourLineHeight.get(contourNeiSegs) is None:
                        s.contourIntersectingEdges.pop()
                        print("Warining! No contour line height defined between: ", contourNeiSegs)
                        continue
                    contourHeight = s.contourLineHeight[
                        s.segmentationData["SegmentationId"][edge[0]], s.segmentationData["SegmentationId"][
                            edge[1]]]
                    contoutUpDownSegments = (
                    s.segmentationData["SegmentationId"][edge[0]], s.segmentationData["SegmentationId"][edge[1]])
                else:
                    contourNeiSegs = (
                    s.segmentationData["SegmentationId"][edge[1]], s.segmentationData["SegmentationId"][edge[0]])
                    if s.contourLineHeight.get(contourNeiSegs) is None:
                        s.contourIntersectingEdges.pop()
                        print("Warining! No contour line height defined between: ", contourNeiSegs)

                        continue
                    contourHeight = s.contourLineHeight[
                        s.segmentationData["SegmentationId"][edge[1]], s.segmentationData["SegmentationId"][
                            edge[0]]]
                    contoutUpDownSegments = (
                    s.segmentationData["SegmentationId"][edge[1]], s.segmentationData["SegmentationId"][edge[0]])

                # check if h1 == h2
                if h1 == h2:
                    w1 = 0.5
                else:
                    w1 = (contourHeight - h2) / (h1 - h2)

                assert 0 <= w1 <= 1
                w2 = 1 - w1

                s.contourLineConstraintWeight.append((w1, w2))
                s.contourLineConstraintHeight.append(contourHeight)
                s.contourLineConstraintUpDownSegments.append(contoutUpDownSegments)

            if len(s.contourIntersectingEdges) != len(s.contourLineConstraintWeight):
                break

    def getUnorderedContourLineConstraintsAtHeight(s, saddleId):
        assert s.nodes[saddleId].criticalType == s.saddleTypeId

        scalarField2D = s.segmentationData[s.segmentationDataScalarName].reshape(s.gridSize)

        allContourEdges = []
        directionsForEdges = [
            (1, 0),
            (0, 1),
        ]

        height = s.nodes[saddleId].scalar

        for d in directionsForEdges:
            s1 = scalarField2D[:scalarField2D.shape[0] - d[0], :scalarField2D.shape[1] - d[1]]
            s2 = scalarField2D[d[0]:, d[1]:]

            contourT1 = np.logical_and(s1 >= height, s2 <= height)
            contourT2 = np.logical_and(s1 <= height, s2 >= height)

            contourAll = np.logical_or(contourT1, contourT2)

            contouredges2D = np.where(contourAll)
            contourEdgesV1Flatten = flatten2DIndex(contouredges2D[0], contouredges2D[1], s.gridSize)
            contourEdgesV2Flatten = flatten2DIndex(contouredges2D[0] + d[0], contouredges2D[1] + d[1], s.gridSize)

            allContourEdges.extend(
                [(contourEdgesV1Flatten[iE], contourEdgesV2Flatten[iE]) for iE in
                 range(contourEdgesV1Flatten.shape[0])])

        contourLineConstraintWeight = []
        contourLineConstraintHeight = []
        contourLineAllPts = []

        # extract the edges that went through all the height
        for iContourEdge in range(len(allContourEdges)):
            edge = allContourEdges[iContourEdge]

            h1 = s.segmentationData[s.segmentationDataScalarName][edge[0]]
            h2 = s.segmentationData[s.segmentationDataScalarName][edge[1]]

            contourLineConstraintHeight.append(height)

            # check if h1 == h2
            if h1 == h2:
                w1 = 0.5
            else:
                w1 = (height - h2) / (h1 - h2)

            w2 = 1 - w1

            contourLineConstraintWeight.append((w1, w2))

            pts = np.array(to2DIndex(edge[0], s.gridSize)) * w1 + np.array(to2DIndex(edge[1], s.gridSize)) * w2
            contourLineAllPts.append(pts)

        # plot
        # plt.figure()
        # contourLineAllPts = np.array(contourLineAllPts)
        # plt.scatter(contourLineAllPts[:,0], contourLineAllPts[:,1])
        # plt.show()

        return allContourEdges, contourLineConstraintWeight, contourLineConstraintHeight

    def extractContourLineConstraints3(s, draw=False,  waitTime = 0):
        '''
        extract the contourline directly from the saddle position
        avoid using the segmentation because they may not form a correct singly connected component
        '''

        for iNode in range(len(s.nodes)):
            # skip the nodes that are not saddle points
            if s.nodes[iNode].criticalType != s.saddleTypeId:
                continue

            allContourEdges, contourLineConstraintWeight, contourLineConstraintHeight = s.getUnorderedContourLineConstraintsAtHeight(
                iNode)

            print("Saddle point id: ", iNode)

            iMeshVert = s.nodeToField[iNode]
            edgesToRemove = []
            saddleNeighborEdges = findSaddleNeighborhood(iMeshVert, s.gridSize)
            initPEdge, initCEdge, currentEdgeIndex = findStartingEdges(iMeshVert, saddleNeighborEdges, s.gridSize,
                                                                       allContourEdges)
            newContourConstraints = ContourConstraint(s.gridSize)

            while initCEdge is not None:
                edgesToRemove.append(initCEdge)
                newEdges, newWeights, newHeights = reorderContourPointsOneLoop(saddleNeighborEdges, initPEdge,
                                                                               initCEdge, currentEdgeIndex, s.gridSize,
                                                                               allContourEdges,
                                                                               contourLineConstraintWeight,
                                                                               contourLineConstraintHeight)
                edgesToRemove.append(newEdges[-1])
                initPEdge, initCEdge, currentEdgeIndex = findStartingEdges(iMeshVert, saddleNeighborEdges, s.gridSize,
                                                                           allContourEdges, edgesToRemove)

                newContour = ContourLine(s.gridSize)
                newContour.saddleAllContourEdges = newEdges
                newContour.saddleAllContourWeights = newWeights
                newContour.saddleAllContourHeights = newHeights

                # check if ccw
                if not newContour.checkCCW():
                    print("Orientation is Not CCW, reverse it!")
                    newContour.reverseOrientation()

                newContourConstraints.addContour(newContour)

            # determine the which upper node belongs to with contour
            upperNodes = []

            for iContour in range(newContourConstraints.numContours()):
                newContourConstraints.getContour(iContour).intialize()
                contourLine = newContourConstraints.getContour(iContour)
                contourLine.parameterize()

                for upperNodeId in s.nodes[iNode].upNodes:
                    upperNodeCorrespondingMeshPt = s.nodeToField[upperNodeId]
                    upperNodePos = np.array(to2DIndex(upperNodeCorrespondingMeshPt, s.gridSize))
                    upperNodePosSp = Point(upperNodePos[0], upperNodePos[1])

                    # if newCountour.includePoint(iContour, upperNodePos):
                    if contourLine.shapelyGeometry.contains(upperNodePosSp):
                        contourLine.embracingHigherNodeId = upperNodeId
                        upperNodes.append(upperNodePos)

            assert newContourConstraints.numContours() == 2

            print("Num countour lines for saddle ", iNode, ":", newContourConstraints.numContours())
            s.saddleContours[iNode] = newContourConstraints

            # determining which contourline contains which higher nodes

            if draw:
                plotSaddleCountourLine(newContourConstraints, s.gridSize, upperNodes=np.array(upperNodes))
                plt.waitforbuttonpress(waitTime)



    def reOrderContourline(s, draw=False, waitTime = 0):

        for iNode in range(len(s.nodes)):
            # skip the nodes that are not saddle points
            if s.nodes[iNode].criticalType != s.saddleTypeId:
                continue
            print("Saddle point id: ", iNode)

            iMeshVert = s.nodeToField[iNode]
            edgesToRemove = []
            saddleNeighborEdges = findSaddleNeighborhood(iMeshVert, s.gridSize)
            initPEdge, initCEdge, currentEdgeIndex = findStartingEdges(iMeshVert, saddleNeighborEdges, s.gridSize,
                                                                       s.contourIntersectingEdges)
            newContourConstraints = ContourConstraint(s.gridSize)


            while initCEdge is not None:
                edgesToRemove.append(initCEdge)
                newEdges, newWeights, newHeights = reorderContourPointsOneLoop(saddleNeighborEdges, initPEdge,
                                                                               initCEdge, currentEdgeIndex, s.gridSize,
                                                                               s.contourIntersectingEdges, s.contourLineConstraintWeight,
                                                                               s.contourLineConstraintHeight)
                edgesToRemove.append(newEdges[-1])
                initPEdge, initCEdge, currentEdgeIndex = findStartingEdges(iMeshVert, saddleNeighborEdges, s.gridSize,
                                                                           s.contourIntersectingEdges, edgesToRemove)


                newContour = ContourLine(s.gridSize)
                newContour.saddleAllContourEdges = newEdges
                newContour.saddleAllContourWeights = newWeights
                newContour.saddleAllContourHeights = newHeights

                # check if ccw
                if not newContour.checkCCW():
                    print("Orientation is Not CCW, reverse it!")
                    newContour.reverseOrientation()

                newContourConstraints.addContour(newContour)

            # determine the which upper node belongs to with contour
            upperNodes = []

            for iContour in range(newContourConstraints.numContours()):
                newContourConstraints.getContour(iContour).intialize()
                contourLine = newContourConstraints.getContour(iContour)
                contourLine.parameterize()

                for upperNodeId in s.nodes[iNode].upNodes:
                    upperNodeCorrespondingMeshPt = s.nodeToField[upperNodeId]
                    upperNodePos = np.array(to2DIndex(upperNodeCorrespondingMeshPt, s.gridSize))
                    upperNodePosSp = Point(upperNodePos[0], upperNodePos[1])

                    # if newCountour.includePoint(iContour, upperNodePos):
                    if contourLine.shapelyGeometry.contains(upperNodePosSp):
                        contourLine.embracingHigherNodeId = upperNodeId
                        upperNodes.append(upperNodePos)

            assert newContourConstraints.numContours() == 2

            print("Num countour lines for saddle ", iNode, ":", newContourConstraints.numContours())
            s.saddleContours[iNode] = newContourConstraints

            # determining which contourline contains which higher nodes


            if draw:
                plotSaddleCountourLine(newContourConstraints, s.gridSize, upperNodes=np.array(upperNodes))
                plt.waitforbuttonpress(waitTime)
    def findRootNode(s):
        for iNode in range(len(s.nodes)):
            if len(s.node(iNode).downNodes) == 0:
                s.rootNodeId = iNode
                break
        assert s.rootNodeId is not None
        return




def findSquareFromEdge(previousEdge, currentEdge, gridSize):
    """
    Compute the two squares adjacent to the current edge and determine the square to use for this step.
    previousEdge: [vertex0, vertex1] (flattened indices), the edge considered from the previous step
    currentEdge: [vertex0, vertex1] (flattened indices), the edge considered in the current step
    Output: three edges in the square excluding the current edge
    """

    i_0, j_0 = to2DIndex(currentEdge[0], gridSize)
    i_1, j_1 = to2DIndex(currentEdge[1], gridSize)

    assert (i_0 - i_1 == 0 and abs(j_0 - j_1) == 1) or (abs(i_0 - i_1) == 1 and j_0 - j_1 == 0)
    # horizontal edge of a square
    if i_0 - i_1 == 0 and abs(j_0 - j_1) == 1:
        square_0 = [(flatten2DIndex(i_0, j_0, gridSize), flatten2DIndex(i_0+1, j_0, gridSize)),
                    (flatten2DIndex(i_0+1, j_0, gridSize), flatten2DIndex(i_1+1, j_1, gridSize)),
                    (flatten2DIndex(i_1+1, j_1, gridSize), flatten2DIndex(i_1, j_1, gridSize))]
        square_1 = [(flatten2DIndex(i_0, j_0, gridSize), flatten2DIndex(i_0 - 1, j_0, gridSize)),
                    (flatten2DIndex(i_0 - 1, j_0, gridSize), flatten2DIndex(i_1 - 1, j_1, gridSize)),
                    (flatten2DIndex(i_1 - 1, j_1, gridSize), flatten2DIndex(i_1, j_1, gridSize))]
    # vertical edge of a square
    elif abs(i_0 - i_1) == 1 and j_0 - j_1 == 0:
        square_0 = [(flatten2DIndex(i_0, j_0, gridSize), flatten2DIndex(i_0, j_0-1, gridSize)),
                    (flatten2DIndex(i_0, j_0-1, gridSize), flatten2DIndex(i_1, j_1-1, gridSize)),
                    (flatten2DIndex(i_1, j_1-1, gridSize), flatten2DIndex(i_1, j_1, gridSize))]
        square_1 = [(flatten2DIndex(i_0, j_0, gridSize), flatten2DIndex(i_0, j_0 + 1, gridSize)),
                    (flatten2DIndex(i_0, j_0 + 1, gridSize), flatten2DIndex(i_1, j_1 + 1, gridSize)),
                    (flatten2DIndex(i_1, j_1 + 1, gridSize), flatten2DIndex(i_1, j_1, gridSize))]

    if previousEdge in square_0 or (previousEdge[1], previousEdge[0]) in square_0:
        return square_1
    elif previousEdge in square_1 or (previousEdge[1], previousEdge[0]) in square_1:
        return square_0
    else:
        raise ValueError("The previous edge is not in either square that contains the current edge.")


def checkCurrentSquareIntersections(currentSquare, contourEdges):
    intersections = [edge for edge in currentSquare if edge in contourEdges]
    reverseSquare = [(edge[1], edge[0]) for edge in currentSquare]
    reverseIntersections = [edge for edge in reverseSquare if edge in contourEdges]
    if len(intersections) + len(reverseIntersections) > 1:
        raise ValueError("Too many intersection points in the current square.")
    elif len(intersections) + len(reverseIntersections) == 0:
        raise ValueError("No intersection point in the current square.")
    else:
        return True


def isEdge(edge, gridSize):
    i_0, j_0 = to2DIndex(edge[0], gridSize)
    i_1, j_1 = to2DIndex(edge[1], gridSize)

    if check2DCoordValidility(i_0, j_0, gridSize) and check2DCoordValidility(i_1, j_1, gridSize):
        if (i_0 - i_1 == 0 and abs(j_0 - j_1) == 1) or (abs(i_0 - i_1) == 1 and j_0 - j_1 == 0):
            return True
        else:
            return False
    else:
        raise ValueError("At least one of the vertices is not valid.")


def findSaddleNeighborhood(saddleVert, gridSize):
    """
    Find the neighboring edges of a saddle point
    Output: list of 12 edges around a saddle
    """
    i, j = to2DIndex(saddleVert, gridSize)

    left = flatten2DIndex(i, j - 1, gridSize)
    topleft = flatten2DIndex(i + 1, j - 1, gridSize)
    top = flatten2DIndex(i + 1, j, gridSize)
    topright = flatten2DIndex(i + 1, j + 1, gridSize)
    right = flatten2DIndex(i, j + 1, gridSize)
    bottomright = flatten2DIndex(i - 1, j + 1, gridSize)
    bottom = flatten2DIndex(i - 1, j, gridSize)
    bottomleft = flatten2DIndex(i - 1, j - 1, gridSize)

    saddleNeighborEdges = [(left, topleft), (topleft, top), (top, topright), (topright, right), (right, bottomright),
                           (bottomright, bottom), (bottom, bottomleft), (bottomleft, left)]

    return saddleNeighborEdges


def findStartingEdges(saddleVert, saddleNeighborEdges, gridSize, contourEdges, edgesToRemove=None):
    """
    Find the starting currentEdge and starting previousEdge for one saddle
    saddlePoint: a vertex (int) on the grid
    saddleNeighborEdges: 8 edges in the neighborhood of a saddle point but not incident to the saddle point
    edgesToRemove: the edges in the neighborhood that are already starting or ending of the a loop, remove these
    Output: currentEdge - the current edge to start with, previousEdge - the edge used to initialize the process
    """

    if edgesToRemove is not None:
        saddleNeighborEdges = [edge for edge in saddleNeighborEdges if edge not in edgesToRemove]
        reverseEdgesToRemove = [(edge[1], edge[0]) for edge in edgesToRemove]
        saddleNeighborEdges = [edge for edge in saddleNeighborEdges if edge not in reverseEdgesToRemove]

    # determine the first edge to start with
    for iSaddleEdge in saddleNeighborEdges:
        currentEdge = None
        iSaddleEdge_0 = to2DIndex(iSaddleEdge[0], gridSize)
        iSaddleEdge_1 = to2DIndex(iSaddleEdge[1], gridSize)
        if check2DCoordValidility(iSaddleEdge_0[0], iSaddleEdge_0[1], gridSize) and check2DCoordValidility(iSaddleEdge_1[0],iSaddleEdge_1[1], gridSize):
            if iSaddleEdge in contourEdges:
                currentEdge = iSaddleEdge
                break
            elif (iSaddleEdge[1], iSaddleEdge[0]) in contourEdges:
                currentEdge = (iSaddleEdge[1], iSaddleEdge[0])
                break

    print("initial current edge: ", currentEdge)
    # if there are no more intersections in the saddle's neighborhood, move on to the next saddle
    if currentEdge is None:
        return None, None, None

    currentEdgeIndex = contourEdges.index(currentEdge)

    if isEdge((saddleVert, currentEdge[0]), gridSize):
        previousEdge = (saddleVert, currentEdge[0])
    elif isEdge((saddleVert, currentEdge[1]), gridSize):
        previousEdge = (saddleVert, currentEdge[1])
    else:
        raise ValueError("The starting edge is not in the neighborhood of the saddle point.")

    return previousEdge, currentEdge, currentEdgeIndex


def plotEdges(edgesToPlot, step=None):
    """
    Plot the edges up to a specified step.
    """
    x = []
    y = []

    for i, edge in enumerate(edgesToPlot):
        if step is not None and i == step:
            break
        i_0, j_0 = to2DIndex(edge[0], gridSize)
        i_1, j_1 = to2DIndex(edge[1], gridSize)
        x.append(j_0)
        x.append(j_1)
        y.append(i_0)
        y.append(i_1)
    plt.scatter(x, y)
    plt.show()


def reorderContourPointsOneLoop(saddleNeighborEdges, previousEdge, currentEdge, currentEdgeIndex, gridSize, contourEdges, contourWeights, contourHeights):
    """
    Reorder the contour line points, weights and heights for one loop around a saddle point (a saddle can have multiple loops).
    Output: reordered contourEdges, contourWeights, and contourHeights
    """

    if previousEdge == currentEdge or (previousEdge[1], previousEdge[0]) == currentEdge:
        raise ValueError("Previous edge and current edge are the same edge.")

    contourEdgesReordered = [currentEdge]
    contourWeightsReordered = [contourWeights[currentEdgeIndex]]
    contourHeightsReordered = [contourHeights[currentEdgeIndex]]

    # plotEdges([currentEdge, previousEdge])
    currentSquare = findSquareFromEdge(previousEdge, currentEdge, gridSize)

    contourEdgesSet = set([(e[0], e[1]) for e in contourEdges])

    # print(checkCurrentSquareIntersections(currentSquare, contourEdges))
    if checkCurrentSquareIntersections(currentSquare, contourEdges):
        for edge in currentSquare:
            if edge in contourEdges:
                previousEdge = currentEdge
                currentEdge = edge
                currentEdgeIndex = contourEdges.index(currentEdge)
            elif (edge[1], edge[0]) in contourEdges:
                previousEdge = currentEdge
                currentEdge = (edge[1], edge[0])
                currentEdgeIndex = contourEdges.index(currentEdge)
    contourEdgesReordered.append(currentEdge)
    contourWeightsReordered.append(contourWeights[currentEdgeIndex])
    contourHeightsReordered.append(contourHeights[currentEdgeIndex])


    while currentEdge not in saddleNeighborEdges and (currentEdge[1], currentEdge[0]) not in saddleNeighborEdges:
        # print("in while loop.")
        currentSquare = findSquareFromEdge(previousEdge, currentEdge, gridSize)
        # print(checkCurrentSquareIntersections(currentSquare, contourEdges))
        if checkCurrentSquareIntersections(currentSquare, contourEdges):
            for edge in currentSquare:

                if edge in contourEdgesSet:
                    previousEdge = currentEdge
                    currentEdge = edge
                    currentEdgeIndex = contourEdges.index(currentEdge)
                elif (edge[1], edge[0]) in contourEdgesSet:
                    previousEdge = currentEdge
                    currentEdge = (edge[1], edge[0])
                    currentEdgeIndex = contourEdges.index(currentEdge)
        contourEdgesReordered.append(currentEdge)
        contourWeightsReordered.append(contourWeights[currentEdgeIndex])
        contourHeightsReordered.append(contourHeights[currentEdgeIndex])

    return contourEdgesReordered, contourWeightsReordered, contourHeightsReordered

def plotSaddleCountourLine(newContourConstraints, gridSize, upperNodes=None):
    # plt.figure()
    # plt.ylim(0, gridSize[1])
    # plt.xlim(0, gridSize[0])
    assert newContourConstraints.numContours() == 2

    colors = ['r', 'g', ]
    cmaps = ['viridis', 'jet', ]
    fig, ax = plt.subplots()

    ax.set_xlim([0, gridSize[0]])
    ax.set_ylim([0, gridSize[1]])

    for iContour in range(newContourConstraints.numContours()):
        allPts = newContourConstraints.getAllNodes(iContour)

        allPts = np.vstack([allPts, allPts[:1, :]])

        points = np.array([allPts[:,0], allPts[:,1],]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        cols = newContourConstraints.getContour(iContour).contourLineParameters

        if upperNodes is not None:
            ax.scatter(upperNodes[iContour, 0], upperNodes[iContour,1], marker ="*", color = colors[iContour])
        # plt.plot(allPts[:, 0], allPts[:, 1], color = colors[iContour])


        lc = LineCollection(segments, cmap=cmaps[iContour])
        lc.set_array(cols)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)
        fig.colorbar(line, ax=ax)

    plt.show()

def writeOBj(outObj, X, Y, Z, gridSize, ):
    '''
    Write a scalar field to a grid mesh
    '''
    Nx = gridSize[0]
    Ny = gridSize[1]
    file  = open(outObj, 'w')
    for i, j in itertools.product(range(Nx), range(Ny)):
        file.write('v %f %f %f\n' %( X[i,j],  Y[i,j], Z[i,j] ))

    # add the faces to the obj file
    for i, j in itertools.product(range(0, Nx-1), range(1, Ny)):
        vId = j + i * Ny
        file.write('f %d %d %d\n' %(vId, vId+1,  vId+Ny+1, ))
        file.write('f %d %d %d\n' %(vId, vId+Ny+1,  vId+Ny, ))

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
        # print("Min distance: ", dis[matchedVIds[-1]])
        assert dis[matchedVIds[-1]] < 1e-6

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
                    contourEdges.pop()
                    continue
                contourHeight = contourLineHeight[fieldSeg["SegmentationId"][edge[0]], fieldSeg["SegmentationId"][edge[1]]]
            else:
                if contourLineHeight.get(
                        (fieldSeg["SegmentationId"][edge[1]], fieldSeg["SegmentationId"][edge[0]])) is None:
                    contourEdges.pop()
                    continue
                contourHeight = contourLineHeight[
                    fieldSeg["SegmentationId"][edge[1]], fieldSeg["SegmentationId"][edge[0]]]

            # check if h1 == h2
            if h1 == h2:
                w1 = 0.5
            else:
                w1 = (contourHeight - h2) / (h1 - h2)

            assert 0 <= w1 <=1
            w2 = 1 - w1
            contourLineConstraintWeight.append([w1, w2])
            contourLineConstraintHeight.append(contourHeight)
        if len(contourEdges) != len(contourLineConstraintWeight):
            break

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

def contourLineInterpolationConstraints(O, r, gridSize):

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

    objFiles = glob.glob(inObjFolder + r'/*.' + inFileExt)
    print(inObjFolder)
    print(outVtkFolder)

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

        mesh.point_data['Height'] = mesh.points[:, 2]

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


def generate_directions(n):
    """
    Generate an [-n,n] x [-n,n] patch of directions to find the segmentations in the neighborhood of a saddle point.
    output: [(-1, -1), (0, -1), ...] (without (0,0))
    If there is an assertion error from the findContourLineHeight() function, try a larger patch.
    """
    x_direction, y_direction = np.arange(-n, n+1, dtype=int), np.arange(-n, n+1, dtype=int)
    directions = list(itertools.product(x_direction, y_direction))
    directions.remove((0, 0))
    return directions

def findContourLineHeightBasedOnEdge(nodes, edges, fieldSeg, gridSize, ):
    """
    In the merge tree, nodes corresponds to critical points and the edges corresponds to the segments
    """
    contourLineHeight = {} # key (a, b): a: higher segment, c lower segment
    iMeshVerts = matchTreeToGrid(nodes.points, fieldSeg.points)

    for iNode in range(nodes.points.shape[0]):
        # non-saddle node
        if nodes['CriticalType'][iNode] != 2:
            continue
        # saddle node
        # iterate all edges to find two nodes above the node
        upSegs = []
        downsegs = []
        for iEdge in range(edges["upNodeId"].shape[0]):
            if edges["downNodeId"][iEdge] == iNode:
                upSegs.append(edges["SegmentationId"][iEdge])
            if edges["upNodeId"][iEdge] == iNode:
                downsegs.append(edges["SegmentationId"][iEdge])

        assert len(upSegs) == 2
        assert len(downsegs) == 1

        contourLineHeight[(upSegs[0], downsegs[0])] = nodes['Scalar'][iNode]
        contourLineHeight[(upSegs[1], downsegs[0])] = nodes['Scalar'][iNode]

    return contourLineHeight

def findContourLineHeight(nodes, fieldSeg, gridSize, directions):
    contourLineHeight = {} # key (a, b): a: higher segment, c lower segment
    iMeshVerts = matchTreeToGrid(nodes.points, fieldSeg.points)

    for iNode in range(nodes.points.shape[0]):

        if nodes['CriticalType'][iNode] != 2:
            continue
        # how to find the height of a contour line?
        # the height of a contour line is the height of a saddle point where two contour line meets
        # at the neighborhood of each saddle point, there will be three types of nodes: two higher segment and one lower segment
        # we will need to find the two higher seg, say (a,b) and the lower seg c
        # the height of contour line a-c, b-c will be the height (scalar) of this saddle point

        # match node to grid mesh
        iMeshVert = iMeshVerts[iNode]
        i, j = to2DIndex(iMeshVert, gridSize)

        # segsInNeighborhood = []
        segsInNeighborhood_dict = {}
        heights_dict = {}

        for d in directions:
            if check2DCoordValidility(i+d[0], j+d[1], gridSize, ):
                neiVertId = flatten2DIndex(i+d[0], j+d[1], gridSize)
                if fieldSeg['SegmentationId'][neiVertId] not in segsInNeighborhood_dict.keys():
                    # segsInNeighborhood.append(fieldSeg['SegmentationId'][neiVertId])
                    segsInNeighborhood_dict[fieldSeg['SegmentationId'][neiVertId]] = 1
                    heights_dict[fieldSeg['SegmentationId'][neiVertId]] = fieldSeg['Height'][neiVertId]
                    # heights.append(fieldSeg['Height'][neiVertId])
                else:
                    segsInNeighborhood_dict[fieldSeg['SegmentationId'][neiVertId]] += 1
        # print(segsInNeighborhood_dict)
        segsInNeighborhood_dict = dict(sorted(segsInNeighborhood_dict.items(), key=lambda item: item[1]))
        # print(segsInNeighborhood_dict)
        # print(heights_dict)
        segsInNeighborhood = list(segsInNeighborhood_dict.keys())
        if len(segsInNeighborhood) > 3:
            assert segsInNeighborhood_dict[segsInNeighborhood[0]] < segsInNeighborhood_dict[segsInNeighborhood[1]]
        segsInNeighborhood = segsInNeighborhood[-3:]
        heights = [heights_dict[key] for key, _ in segsInNeighborhood_dict.items()]
        heights = heights[-3:]
        # print(heights)
        # print(segsInNeighborhood)
        # print(iNode)
        assert len(segsInNeighborhood) == 3
        sortedId = np.argsort(heights)
        contourLineHeight[(segsInNeighborhood[sortedId[2]], segsInNeighborhood[sortedId[0]])] = nodes['Scalar'][iNode]
        contourLineHeight[(segsInNeighborhood[sortedId[1]], segsInNeighborhood[sortedId[0]])] = nodes['Scalar'][iNode]

        # print(contourLineHeight)

    return contourLineHeight


if __name__ == '__main__':
    inputScalarField = './Data/S04_GenerateNewScalarField.py/M13338.obj'
    inputSegmentedField = './Data/S04_GenerateNewScalarField.py/Topology/M13338/Seg.vtk'
    inputMergeTreeNodesField = './Data/S04_GenerateNewScalarField.py/Topology/M13338/Node.vtk'
    inputEdges = './Data/S04_GenerateNewScalarField.py/Topology/M13338/Edges.vtk'
    gridSize = (200, 200)
    dType = np.float64
    directions = [
        (-1, -1),
        (0, -1),
        (1, -1),
        (-1, 0),
        (1, 0),
        (-1, 1),
        (0, 1),
        (1, 1),
    ]

    numberOfVariable = gridSize[0] * gridSize[1]

    fieldSeg = pv.read(inputSegmentedField)
    fieldSeg.points[:, 2] = fieldSeg['Height']
    field = pv.PolyData(inputScalarField)

    nodes = pv.read(inputMergeTreeNodesField)
    edges = pv.read(inputEdges)

    contourLineHeight = findContourLineHeight(nodes, edges, fieldSeg, gridSize, directions)
    print(contourLineHeight)

