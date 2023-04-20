import numpy as np

from M01_TopologicalExtraction import *
from copy import deepcopy
from matplotlib.collections import LineCollection


def findCorr(iNode, iTree, treeToTreeCorrespondence):
    anotherTreeId = (iTree + 1) % 2

    i = np.where(treeToTreeCorrespondence[:, iTree] == iNode)[0][0]
    return treeToTreeCorrespondence[i, anotherTreeId]



class LinearAnimation:
    """

    """
    def __init__(s, tree0, tree1, correspondences, intermediateTreeHeights):
        s.tree0 = tree0
        s.tree1 = tree1
        s.correspondences = correspondences
        s.intermediateTreeHeights = intermediateTreeHeights

        s.gridSize = tree0.gridSize
        s.numRegistrationCandidates = 100

        s.intermediateTree = Tree()
        assert  tree0.saddleTypeId == tree1.saddleTypeId
        s.intermediateTree.saddleTypeId = tree0.saddleTypeId
        # s.persistence = persistence

        # convert node correspondences
        s.corrsTree0ToTree1 = np.array([-1 for iNode in range(tree0.numNodes())])
        s.corrsTree1ToTree0 = np.array([-1 for iNode in range(tree1.numNodes())])

        for corr in s.correspondences:
            if corr[1] != -1:
                s.corrsTree1ToTree0[corr[1]] = corr[0]
            if corr[0] != -1:
                s.corrsTree0ToTree1[corr[0]] = corr[1]

        s.numSharedVerts = np.where(s.corrsTree0ToTree1 != -1)[0].shape[0]
        s.numUniqueVertsTree0 = np.where(s.corrsTree0ToTree1 == -1)[0].shape[0]
        s.numUniqueVertsTree1 = np.where(s.corrsTree1ToTree0 == -1)[0].shape[0]

        s.tree0ToIntermediateTree = np.array([-1 for iNode in range(tree0.numNodes())])
        s.tree1ToIntermediateTree = np.array([-1 for iNode in range(tree1.numNodes())])

        # Construct the intermediate tree
        # totally
        s.numNodesInterTree = s.numSharedVerts + s.numUniqueVertsTree0 + s.numUniqueVertsTree1

        # generate nodes
        nodes = []

        for iTree0Node in range(tree0.numNodes()):
            node = deepcopy(tree0.node(iTree0Node))

            node.tree0Corr = iTree0Node
            node.tree1Corr = s.corrsTree0ToTree1[iTree0Node]

            s.tree0ToIntermediateTree[iTree0Node] = iTree0Node

            nodes.append(node)

        # find correspondence to t1
        for iTree1Node in range(tree1.numNodes()):

            if s.corrsTree1ToTree0[iTree1Node,] != -1:
                node = nodes[s.corrsTree1ToTree0[iTree1Node]]
                node.tree1Corr = iTree1Node

            else:
                node = TreeNode()

                node.position = tree1.nodes[iTree1Node].position
                node.criticalType = tree1.nodes[iTree1Node].criticalType
                node.scalar = tree1.nodes[iTree1Node].scalar
                node.posInField = tree1.nodes[iTree1Node].posInField  # position in terms of row, col at the scalar field

                node.id = iTree1Node + tree0.numNodes()
                node.tree0Corr = s.corrsTree1ToTree0[iTree1Node]
                node.tree1Corr = iTree1Node
                nodes.append(node)

            s.tree1ToIntermediateTree[iTree1Node] = node.id

        edges = []

        for iTree0Edge in range(tree0.numEdges()):
            edges.append(deepcopy(tree0.edges[iTree0Edge]))

        # tree1 edges
        for iTree1Edge in range(tree1.numEdges()):
            edge = Edge()

            # newEdge.id = iEdge
            # newEdge.segmentId = edgesData["SegmentationId"][iEdge]
            # newEdge.nodes = [edgesData[actualUp + "NodeId"][iEdge], edgesData[actualDown + "NodeId"][iEdge], ]
            # newEdge.upNode = edgesData[actualUp + "NodeId"][iEdge]
            # newEdge.downNode = edgesData[actualDown + "NodeId"][iEdge]
            # s.edges.append(newEdge)
            # # initialize the connectivity infos for up node
            #
            # upNode = s.nodes[newEdge.upNode]
            # upNode.downEdges.append(iEdge)
            # upNode.downNodes.append(edgesData[actualDown + "NodeId"][iEdge])
            #
            # downNode = s.nodes[newEdge.downNode]
            # downNode.upNodes.append(edgesData[actualUp + "NodeId"][iEdge])
            # downNode.upEdges.append(iEdge)

            tree1Edge = tree1.edges[iTree1Edge]

            edge.id = len(edges)
            edge.segmentId = -1
            edge.segmentIdTree1 = tree1Edge.segmentId
            # edge.segmentIdTree0 = # upNode's tree0Corr's down edge
            upNodeTree1 = tree1Edge.upNode

            edge.upNode = s.tree1ToIntermediateTree[upNodeTree1]
            edge.downNode = s.tree1ToIntermediateTree[tree1Edge.downNode]
            edge.nodes = [edge.upNode, edge.downNode]

            upNode = nodes[edge.upNode]
            upNode.downEdges.append(edge.id)
            upNode.downNodes.append(edge.downNode)

            downNode = nodes[edge.downNode]
            downNode.upNodes.append(edge.upNode)
            downNode.upEdges.append(edge.id)

            edges.append(edge)

        s.intermediateTree.initFrom(nodes, edges)

        # register the contour line

        for iNode in range(s.intermediateTree.numNodes()):
            if s.intermediateTree.node(iNode).tree0Corr == -1:
                s.intermediateTree.node(iNode).type = "emerging"
            elif s.intermediateTree.node(iNode).tree1Corr == -1:
                s.intermediateTree.node(iNode).type = "vanishing"
            else:
                s.intermediateTree.node(iNode).type = "preserving"

        for iNode in range(s.intermediateTree.numNodes()):
            # if s.intermediateTree.node(iNode).tree0Corr == -1:
            #     s.intermediateTree.node(iNode).type = "emerging"
            # elif s.intermediateTree.node(iNode).tree1Corr == -1:
            #     s.intermediateTree.node(iNode).type = "vanishing"
            # else:
            #     s.intermediateTree.node(iNode).type = "preserving"

            # if such a node is a saddle node in either tree0 or tree1 we treat it as a saddle
            if s.getTree0CorrespondingNode(iNode).criticalType == s.intermediateTree.saddleTypeId or\
                s.getTree1CorrespondingNode(iNode).criticalType == s.intermediateTree.saddleTypeId:
                s.intermediateTree.node(iNode).criticalType = s.intermediateTree.saddleTypeId

            if s.intermediateTree.node(iNode).criticalType != s.intermediateTree.saddleTypeId:
                continue

            newCountourConstraints = ContourConstraint(s.gridSize)

            # process contour for saddle saddle
            if s.intermediateTree.node(iNode).tree0Corr != -1 and s.intermediateTree.node(iNode).tree1Corr != -1:
                # preserving saddle
                # interleave the vertices in the curve in this case
                # first step find the matching of the contour line based on the contained up node
                contourConstraintTree0 = s.tree0.saddleContours[s.intermediateTree.node(iNode).tree0Corr]
                contourConstraintTree1 = s.tree1.saddleContours[s.intermediateTree.node(iNode).tree1Corr]

                assert s.getTree0CorrespondingNode(iNode).criticalType == s.intermediateTree.saddleTypeId and \
                        s.getTree1CorrespondingNode(iNode).criticalType == s.intermediateTree.saddleTypeId

                # both contour lines are preserved
                contourMatches =[[0,0], [0,1]]

                nextNodeTree0 = contourConstraintTree0.getContour(contourMatches[0][0]).embracingHigherNodeId
                embracingTree0UpNodeInIntermediateTree = s.tree0ToIntermediateTree[s.getClosestUpNodeTree0(nextNodeTree0)]
                nextNodeTree1 = contourConstraintTree1.getContour(contourMatches[0][1]).embracingHigherNodeId
                embracingTree1UpNodeInIntermediateTree = s.tree0ToIntermediateTree[s.getClosestUpNodeTree1(nextNodeTree1)]

                if embracingTree0UpNodeInIntermediateTree != embracingTree1UpNodeInIntermediateTree:
                    contourMatches = [[0, 1], [1, 0]]

                assert s.tree0ToIntermediateTree[s.getClosestUpNodeTree0(contourConstraintTree0.getContour(contourMatches[0][0]).embracingHigherNodeId)] \
                    == s.tree1ToIntermediateTree[s.getClosestUpNodeTree1(contourConstraintTree1.getContour(contourMatches[0][1]).embracingHigherNodeId)]

                assert s.tree0ToIntermediateTree[s.getClosestUpNodeTree0(contourConstraintTree0.getContour(contourMatches[1][0]).embracingHigherNodeId)] \
                    == s.tree1ToIntermediateTree[s.getClosestUpNodeTree1(contourConstraintTree1.getContour(contourMatches[1][1]).embracingHigherNodeId)]


                for iContour in range(len(contourMatches)):
                    newContour = s.blendContourLines(contourConstraintTree0.getContour(contourMatches[iContour][0]),
                                        contourConstraintTree1.getContour(contourMatches[iContour][1]))
                    newContour.type = "preserving"
                    newContour.startState = contourConstraintTree0.getContour(contourMatches[iContour][0])
                    newContour.endState = contourConstraintTree1.getContour(contourMatches[iContour][1])
                    newCountourConstraints.addContour(newContour)

                s.intermediateTree.saddleContours[iNode] = newCountourConstraints

                pass
            elif s.intermediateTree.node(iNode).tree0Corr != -1:
                saddleTree0 = s.tree0.node( s.intermediateTree.node(iNode).tree0Corr)
                # vanishing saddle
                # Determining whether it's one node vanishing or both node vanishing, based on whether both of its upnodes vanish

                preservingCountourId = -1
                preservingUpNode = -1
                contourLIneConstraints = tree0.saddleContours[saddleTree0.id]

                for iContour, contourLine in enumerate(contourLIneConstraints.contourLines):
                    upNode = contourLine.embracingHigherNodeId

                    if findCorr(upNode, 0, s.correspondences) != -1:
                        preservingCountourId = iContour
                        preservingUpNode = upNode

                # then we have to construct the end state of that contour line, but exacting the level set at height saddleTree0.scalar
                # and includes the upNode

                s.intermediateTree.saddleContours[iNode] = ContourConstraint(s.gridSize)

                contourLineAtHeight = s.getContourLineAtHeight(tree0, s.intermediateTreeHeights[saddleTree0.id, -1])

                # assert contourLineAtHeight[0].includePoint(tree1.node(findCorr(preservingUpNode, 0, s.correspondences)).posInField)

                contourLineEndState = contourLineAtHeight[0]
                newCountourConstraints = ContourConstraint(s.gridSize)

                allPts = contourLineEndState.allNodes

                allPts[:, 0] = 1.4 * allPts[:, 0] / s.gridSize[0] - 0.7
                allPts[:, 1] = 1.4 * allPts[:, 1] / s.gridSize[1] - 0.7
                allPts = np.hstack([allPts, np.zeros((allPts.shape[0], 1))])

                pd = pv.PolyData(allPts)
                pd.save("allPts.ply", binary=False)

                for iContour, contourLine in enumerate(contourLIneConstraints.contourLines):

                    if iContour == preservingCountourId:
                        parameterShifts = np.linspace(0,1, s.numRegistrationCandidates,endpoint=False )
                        bestStartNode = None
                        for iStartNode in range(contourLineEndState.numVertices()):
                            contourLineEndState.parameterize(iStartNode)
                            pts=[]
                            ts = np.linspace(0,1,100)
                            for t in ts:
                                pts.append(contourLineEndState.getPosition(t))
                            pts = np.array(pts)
                            # plt.plot(pts[:,0], pts[:,1])
                            # plt.show()

                            x = pts[:,0]
                            y = pts[:,1]
                            cols =  ts

                            points = np.array([x, y]).T.reshape(-1, 1, 2)
                            segments = np.concatenate([points[:-1], points[1:]], axis=1)

                            fig, ax = plt.subplots()
                            lc = LineCollection(segments, cmap='viridis')
                            lc.set_array(cols)
                            lc.set_linewidth(2)
                            line = ax.add_collection(lc)
                            fig.colorbar(line, ax=ax)

                            plt.show()


                        newContour = s.blendContourLines(contourLine, contourLineEndState)
                        newContour.type = "preserving"
                        newContour.startState = contourLine
                        newContour.endState = contourLineEndState


                    else:
                        newContour = deepcopy(contourLine)
                        newContour.startState = contourLine

                        newContour.type = "vanishing"

                    newCountourConstraints.addContour(newContour)
                s.intermediateTree.saddleContours[iNode] = newCountourConstraints

            elif s.intermediateTree.node(iNode).tree1Corr != -1:
                # emerging saddle
                # in this case we just duplicate the contour line from tree 0
                contourLine = deepcopy(tree1.saddleContours[s.intermediateTree.node(iNode).tree1Corr])
                # match the upper node
                contourLine.embracingHigherNodeId = [s.tree1ToIntermediateTree[contourLine.embracingHigherNodeId[0]],
                                                     s.tree1ToIntermediateTree[contourLine.embracingHigherNodeId[1]],
                                                     ]
                s.intermediateTree.saddleContours[iNode] = contourLine
            else:
                assert False
        s.intermediateTree.findRootNode()

    def blendContourLines(s, contourline0, contourline1):
        newContour = ContourLine(s.gridSize)

        # calculate parameters

        i0 = 0
        i1 = 0

        allParameters = []
        correspondence = [] # 0: from tree0, 1: from tree:1, from both trees
        while i0 < contourline0.numVertices() and i1 < contourline1.numVertices():
            if contourline0.contourLineParameters[i0] == contourline0.contourLineParameters[i0]:
                correspondence.append(2)
                allParameters.append(contourline0.contourLineParameters[i0])
                i0 = i0 +1
                i1 = i1 +1
            elif contourline0.contourLineParameters[i0] == contourline0.contourLineParameters[i0]:
                correspondence.append(0)
                allParameters.append(contourline0.contourLineParameters[i0])
                i0 = i0 +1

            else:
                correspondence.append(1)
                allParameters.append(contourline1.contourLineParameters[i1])
                i1 = i1 +1

        newContour.contourLineParameters = allParameters
        newContour.correspondence = correspondence
        return newContour

    def getTree0CorrespondingNode(s, iNode):
        return s.tree0.node(s.intermediateTree.node(iNode).tree0Corr)

    def getTree1CorrespondingNode(s, iNode):
        return s.tree1.node(s.intermediateTree.node(iNode).tree1Corr)

    def getContourLineAtHeight(s, tree, height):
        scalarField2D = tree.segmentationData[tree.segmentationDataScalarName].reshape(s.gridSize)

        allContourEdges = []
        directionsForEdges = [
            (1, 0),
            (0, 1),
        ]
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
                [(contourEdgesV1Flatten[iE], contourEdgesV2Flatten[iE]) for iE in range(contourEdgesV1Flatten.shape[0])])

        contourLineConstraintWeight = []
        contourLineConstraintHeight = []
        contourLineAllPts = []

        for iContourEdge in range(len(allContourEdges)):
            edge = allContourEdges[iContourEdge]

            h1 =  tree.segmentationData[tree.segmentationDataScalarName][edge[0]]
            h2 =  tree.segmentationData[tree.segmentationDataScalarName][edge[1]]

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

        # Connect Contours
        contourEdgesReordered, contourWeightsReordered, contourHeightsReordered = s.constructNonSaddleContours( allContourEdges, contourLineConstraintWeight, contourLineConstraintHeight)

        newContour = ContourLine(s.gridSize)
        newContour.saddleAllContourEdges = contourEdgesReordered
        newContour.saddleAllContourWeights = contourWeightsReordered
        newContour.saddleAllContourHeights = contourHeightsReordered
        # check if ccw
        if not newContour.checkCCW():
            print("Orientation is Not CCW, reverse it!")
            newContour.reverseOrientation()

        contourLines = [newContour]
        # plt.figure()
        # plt.ylim(0, s.gridSize[1])
        # plt.xlim(0, s.gridSize[0])
        #
        # colors = ['r', 'g', ]
        # allPts = newContour.getAllNodes()
        #
        # allPts = np.vstack([allPts, allPts[:1, :]])
        # plt.plot(allPts[:, 0], allPts[:, 1],)
        #
        # plt.show()

        return contourLines


    def findStartingEdgesForNonSaddleContour(s,currentEdge, allContourEdges):
        """
        Find the starting currentEdge and starting previousEdge for one saddle
        Output: currentEdge - the current edge to start with, previousEdge - the edge used to initialize the process
        """

        # if edgesToRemove is not None:
        #     saddleNeighborEdges = [edge for edge in saddleNeighborEdges if edge not in edgesToRemove]
        #     reverseEdgesToRemove = [(edge[1], edge[0]) for edge in edgesToRemove]
        #     saddleNeighborEdges = [edge for edge in saddleNeighborEdges if edge not in reverseEdgesToRemove]
        #
        # # determine the first edge to start with
        # for iSaddleEdge in saddleNeighborEdges:
        #     currentEdge = None
        #     iSaddleEdge_0 = to2DIndex(iSaddleEdge[0], gridSize)
        #     iSaddleEdge_1 = to2DIndex(iSaddleEdge[1], gridSize)
        #     if check2DCoordValidility(iSaddleEdge_0[0], iSaddleEdge_0[1], gridSize) and check2DCoordValidility(
        #             iSaddleEdge_1[0], iSaddleEdge_1[1], gridSize):
        #         if iSaddleEdge in contourEdges:
        #             currentEdge = iSaddleEdge
        #             break
        #         elif (iSaddleEdge[1], iSaddleEdge[0]) in contourEdges:
        #             currentEdge = (iSaddleEdge[1], iSaddleEdge[0])
        #             break
        #
        # print("initial current edge: ", currentEdge)
        # # if there are no more intersections in the saddle's neighborhood, move on to the next saddle
        # if currentEdge is None:
        #     return None, None, None
        #
        # currentEdgeIndex = contourEdges.index(currentEdge)
        #
        # if isEdge((saddleVert, currentEdge[0]), gridSize):
        #     previousEdge = (saddleVert, currentEdge[0])
        # elif isEdge((saddleVert, currentEdge[1]), gridSize):
        #     previousEdge = (saddleVert, currentEdge[1])
        # else:
        #     raise ValueError("The starting edge is not in the neighborhood of the saddle point.")
        #
        # return previousEdge, currentEdge, currentEdgeIndex

        edgeEnds = (np.array(to2DIndex(currentEdge[0], s.gridSize)),
                    np.array(to2DIndex(currentEdge[1], s.gridSize)))

        if abs(edgeEnds[0][0] - edgeEnds[1][0]) == 1:
            smallerEnd = 0 if edgeEnds[0][0] < edgeEnds[1][0] else 1

            edges1 = s.getSquareEdges(edgeEnds[smallerEnd])
            edges2 = s.getSquareEdges(edgeEnds[smallerEnd] + np.array([0, -1]))

        else:
            smallerEnd = 0 if edgeEnds[0][1] < edgeEnds[1][1] else 1

            edges1 = s.getSquareEdges(edgeEnds[smallerEnd])
            edges2 = s.getSquareEdges(edgeEnds[smallerEnd] + np.array([-1, 0]))

        edgesAllDul = [(flatten2DIndex(e[0][0], e[0][1], s.gridSize), flatten2DIndex(e[1][0], e[1][1], s.gridSize)) for e in edges1] \
                      + [(flatten2DIndex(e[0][0], e[0][1], s.gridSize), flatten2DIndex(e[1][0], e[1][1], s.gridSize)) for e in edges2]

        edgesAll = []

        for e in edgesAllDul:
            if not s.isSameEdge(e, currentEdge):
                if (e[0] > [e[1]]):
                    edgesAll.append((e[1], e[0]))
                else:
                    edgesAll.append(e)
        # then find another edge

        prevEdge = None
        for e in edgesAll:
            if e in allContourEdges or (e[1], e[0]) in allContourEdges:
                prevEdge = e
                break

        assert prevEdge is not None

        return prevEdge


    def getSquareEdges(s, topLeft):
        p = np.array(topLeft)

        d = np.array([(1, 0), (0, 1), (1, 1)])

        edges = [(p, p + d[0]),
                 (p, p + d[1]),
                 (p+ d[0], p + d[2]),
                 (p+ d[1], p + d[2]),
                 ]

        return edges

    def isSameEdge(s, e1, e2):
        """
        e1, e2: edges in flattened index
        """

        if e1[0] == e2[0] and e1[1] == e2[1]:
            return True

        elif e1[0] == e2[1] and e1[1] == e2[0]:
            return True
        return  False

    def constructNonSaddleContours(s, allContourEdges, contourLineConstraintWeight, contourLineConstraintHeight):
        currentEdgeIndex = 0
        currentEdge = allContourEdges[currentEdgeIndex]

        startEdge = currentEdge

        previousEdge = s.findStartingEdgesForNonSaddleContour(currentEdge, allContourEdges)

        currentSquare = findSquareFromEdge(previousEdge, currentEdge, s.gridSize)

        contourEdgesSet = set([(e[0], e[1]) for e in allContourEdges])

        contourEdgesReordered = [currentEdge]
        contourWeightsReordered = [contourLineConstraintWeight[currentEdgeIndex]]
        contourHeightsReordered = [contourLineConstraintHeight[currentEdgeIndex]]

        while len(contourEdgesReordered) != len(allContourEdges):
            # print("in while loop.")
            currentSquare = findSquareFromEdge(previousEdge, currentEdge, s.gridSize)
            # print(checkCurrentSquareIntersections(currentSquare, contourEdges))
            if checkCurrentSquareIntersections(currentSquare, allContourEdges):
                for edge in currentSquare:

                    if edge in contourEdgesSet:
                        previousEdge = currentEdge
                        currentEdge = edge
                        currentEdgeIndex = allContourEdges.index(currentEdge)
                    elif (edge[1], edge[0]) in contourEdgesSet:
                        previousEdge = currentEdge
                        currentEdge = (edge[1], edge[0])
                        currentEdgeIndex = allContourEdges.index(currentEdge)
            contourEdgesReordered.append(currentEdge)
            contourWeightsReordered.append(contourLineConstraintWeight[currentEdgeIndex])
            contourHeightsReordered.append(contourLineConstraintHeight[currentEdgeIndex])

        return contourEdgesReordered, contourWeightsReordered, contourHeightsReordered

    def getTree0NodePersistanceType(s, iNode):

        return s.intermediateTree.node(s.tree0ToIntermediateTree[iNode]).type


    def getTree1NodePersistanceType(s, iNode):
        return s.intermediateTree.node(s.tree1ToIntermediateTree[iNode]).type

    # this is incomplete, assuming the one level higher node is a preserving node
    def getClosestUpNodeTree0(s, nextNodeTree):

        while s.intermediateTree.node(s.tree0ToIntermediateTree[nextNodeTree]).type != "preserving":
            upNodes = s.tree0.node(nextNodeTree).upNodes

            assert s.intermediateTree.node(
                s.tree0ToIntermediateTree[nextNodeTree]).criticalType == s.intermediateTree.saddleTypeId
            if s.getTree0NodePersistanceType(upNodes[0]) == "preserving":
                nextNodeTree = upNodes[0]
            else:
                nextNodeTree = upNodes[1]

        return nextNodeTree

    def getClosestUpNodeTree1(s, nextNodeTree):

        while s.intermediateTree.node(s.tree1ToIntermediateTree[nextNodeTree]).type != "preserving":
            upNodes = s.tree1.node(nextNodeTree).upNodes

            assert s.intermediateTree.node(
                s.tree1ToIntermediateTree[nextNodeTree]).criticalType == s.intermediateTree.saddleTypeId
            if s.getTree1NodePersistanceType(upNodes[0]) == "preserving":
                nextNodeTree = upNodes[0]
            else:
                nextNodeTree = upNodes[1]

        return nextNodeTree


    def interpolateAnimation(s, t):

        nodeQueue = [s.intermediateTree.rootNodeId]

        while len(nodeQueue) != 0:
            nodeId = nodeQueue[0]
            node = s.intermediateTree.node(nodeId)

            if node.type == "preserving":
                # interpolate
                pass

            # elif node.type == ""


    def computeRegisterationDis(s, contourLine1, contourLine2):

        for iVert in range(len(contourLine1.numVertices())):
            # p1 =
            pass
