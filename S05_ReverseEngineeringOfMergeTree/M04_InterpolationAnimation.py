from M01_TopologicalExtraction import *
from copy import deepcopy
class LinearAnimation:
    """

    """
    def __init__(s, tree0, tree1, correspondences):
        s.tree0 = tree0
        s.tree1 = tree1
        s.correspondences = correspondences

        s.gridSize = tree0.gridsize

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


        for iTree1Node in range(tree1.numNodes()):
            node = TreeNode()

            node.position = tree1.nodes[iTree1Node].position
            node.criticalType = tree1.nodes[iTree1Node].criticalType
            node.scalar = tree1.nodes[iTree1Node].scalar
            node.posInField = tree1.nodes[iTree1Node].posInField  # position in terms of row, col at the scalar field

            node.id = iTree1Node + tree0.numNodes()
            node.tree0Corr = s.corrsTree1ToTree0[iTree1Node]
            node.tree1Corr = iTree1Node

            s.tree1ToIntermediateTree[iTree1Node] = node.id
            nodes.append(node)

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

            # if such a node is a saddle node in either tree0 or tree1 we treat it as a saddle
            if s.getTree0CorrespondingNode(iNode).criticalType == s.intermediateTree.saddleTypeId or\
                s.getTree1CorrespondingNode(iNode).criticalType == s.intermediateTree.saddleTypeId:
                s.intermediateTree.node(iNode).criticalType = s.intermediateTree.saddleTypeId

            if s.intermediateTree.node(iNode).criticalType != s.intermediateTree.saddleTypeId:
                continue

            newCountourConstraints = CountourConstraint(gridSize)

            # process contour for saddle saddle
            if s.intermediateTree.node(iNode).tree0Corr != -1 and s.intermediateTree.node(iNode).tree1Corr != -1:
                # preserving saddle
                # interleave the vertices in the curve in this case
                # first step find the matching of the contour line based on the contained up node
                contourConstraintTree0 = s.tree0.saddleContours[s.intermediateTree.node(iNode).tree0Corr]
                contourConstraintTree1 = s.tree1.saddleContours[s.intermediateTree.node(iNode).tree1Corr]

                if s.getTree0CorrespondingNode(iNode).criticalType == s.intermediateTree.saddleTypeId or \
                        s.getTree1CorrespondingNode(iNode).criticalType == s.intermediateTree.saddleTypeId:
                    # both contour lines are preserved
                    contourMatches =[[0,0], [0,1]]
                    if s.tree0ToIntermediateTree[contourConstraintTree0.getContour(0).embracingHigherNodeId] \
                        != s.tree1ToIntermediateTree[contourConstraintTree1.getContour(0).embracingHigherNodeId]:
                        contourMatches = [[0, 1], [1, 0]]

                    assert s.tree0ToIntermediateTree[contourConstraintTree0.getContour(contourMatches[0][0]).embracingHigherNodeId] \
                        == s.tree1ToIntermediateTree[contourConstraintTree1.getContour(contourMatches[0][1]).embracingHigherNodeId]

                    assert s.tree0ToIntermediateTree[contourConstraintTree0.getContour(contourMatches[1][0]).embracingHigherNodeId] \
                        == s.tree1ToIntermediateTree[contourConstraintTree1.getContour(contourMatches[1][1]).embracingHigherNodeId]


                    for iContour in range(len(contourMatches)):
                        newContour = s.blendContourLines(contourConstraintTree0.getContour(contourMatches[iContour][0]),
                                            contourConstraintTree0.getContour(contourMatches[iContour][0]))
                        newContour.type = "preserving"
                        newCountourConstraints.addContour(newContour)
                else:
                    

                pass
            elif s.intermediateTree.node(iNode).tree0Corr != -1:
                # in this case we just duplicate the contour line from tree 0
                contourLine = deepcopy(tree0.saddleContours[s.intermediateTree.node(iNode).tree0Corr])
                # match the upper node
                contourLine.embracingHigherNodeId = [s.tree0ToIntermediateTree[contourLine.embracingHigherNodeId[0]],
                                                     s.tree0ToIntermediateTree[contourLine.embracingHigherNodeId[1]],
                                                     ]
                s.intermediateTree.saddleContours[iNode] = contourLine
                pass
            elif s.intermediateTree.node(iNode).tree1Corr != -1:
                # in this case we just duplicate the contour line from tree 0
                contourLine = deepcopy(tree1.saddleContours[s.intermediateTree.node(iNode).tree1Corr])
                # match the upper node
                contourLine.embracingHigherNodeId = [s.tree1ToIntermediateTree[contourLine.embracingHigherNodeId[0]],
                                                     s.tree1ToIntermediateTree[contourLine.embracingHigherNodeId[1]],
                                                     ]
                s.intermediateTree.saddleContours[iNode] = contourLine
            else:
                assert False

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





