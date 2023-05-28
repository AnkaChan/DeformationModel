from M01_TopologicalExtraction import *

import pyvista as pv
from os.path import join
import vtk
import numpy as np
from qpsolvers import solve_qp, available_solvers
from scipy import sparse

class ContourLineGridIntersector:
    def __init__(s, gridSize):
        s.gridSize = gridSize
        s.equalityConstraintIds = []
        s.equalityConstraintWeights = []
        s.equalityConstraintVal = []

        s.fig, s.ax = plt.subplots()
        s.scatter = s.ax.scatter([0], [0], s=0.1)
        s.ax.set_xlim([0, s.gridSize[0]])
        s.ax.set_ylim([0, s.gridSize[1]])
        s.ax.grid()

    def calculate_intersection(s, seg1, seg2):
        x1, y1 = seg1[0]
        x2, y2 = seg1[1]
        x3, y3 = seg2[0]
        x4, y4 = seg2[1]

        # Calculate the denominator
        denominator = ((x1 - x2) * (y3 - y4)) - ((y1 - y2) * (x3 - x4))

        # Calculate the numerator for x and y coordinates
        num_x = ((x1 * y2 - y1 * x2) * (x3 - x4)) - ((x1 - x2) * (x3 * y4 - y3 * x4))
        num_y = ((x1 * y2 - y1 * x2) * (y3 - y4)) - ((y1 - y2) * (x3 * y4 - y3 * x4))

        # Calculate the intersection point
        if denominator != 0:
            intersection_x = num_x / denominator
            intersection_y = num_y / denominator
            return intersection_x, intersection_y
        else:
            # Lines are parallel or coincident, no intersection
            return None
    def contourLineArea(s, contourLine):
        startingP = np.array(contourLine[0])

        totalOrientedArea = 0

        for i in range(1, len(contourLine)-1):
            edgeP1 = np.array(contourLine[i])
            iEdge = i
            if iEdge == len(contourLine) - 1:
                nextNode = 0
            else:
                nextNode = iEdge + 1

            nextEdge = np.array(contourLine[nextNode]) - edgeP1

            p0p1 = startingP - edgeP1

            totalOrientedArea = totalOrientedArea + nextEdge[0] * p0p1[1] -  nextEdge[1] * p0p1[0]

        assert totalOrientedArea > 0

        return totalOrientedArea

    def intersectWithGrid(s, inputContourLineData, plot=False):
        # "ContourLines"
        # "ContourLineHeights"
        # "Extremity"

        s.equalityConstraintIds.clear()
        s.equalityConstraintWeights.clear()
        s.equalityConstraintVal.clear()

        s.minContourArea = 20

        # extremities
        for iNode in range(len(inputContourLineData["Extremity"])):
            node = inputContourLineData["Extremity"][iNode]
            pos = [int(node[0]), int(node[1])]
            s.equalityConstraintIds.append([int(flatten2DIndex(*pos, s.gridSize))])
            s.equalityConstraintWeights.append([1])
            s.equalityConstraintVal.append(node[2])

        # contour lines
        for iContour in range(len(inputContourLineData["ContourLines"])):
            contourLine = inputContourLineData["ContourLines"][iContour]

            contourLineArea = s.contourLineArea(contourLine)
            if contourLineArea < s.minContourArea:
                print("Skipping contour: ", iContour, " because of its area being too small.")
                continue

            for iEdge in range(len(contourLine)):
                iNode1 = iEdge
                iNode2 = iEdge + 1 if iEdge < len(contourLine)-1 else 0
                lineSeg = np.array([
                    contourLine[iNode1],
                    contourLine[iNode2]
                ])

                xMin = int(np.ceil(min(lineSeg[0,0], lineSeg[1,0])))
                xMax = int(np.floor(max(lineSeg[0,0], lineSeg[1,0])))

                yMin = int(np.ceil(min(lineSeg[0,1], lineSeg[1,1])))
                yMax = int(np.floor(max(lineSeg[0,1], lineSeg[1,1])))

                # if xMax < xMin or yMax < yMin:
                #     continue

                xs = np.linspace(xMin, xMax, num= (1+xMax-xMin), endpoint=True)

                ys = np.linspace(yMin, yMax, num= (1+yMax-yMin), endpoint=True)

                for x in xs:
                    scanline = np.array([
                        [x, yMin-1],
                        [x, yMax+1],
                    ])

                    intersectRes = s.calculate_intersection(scanline, lineSeg)
                    if intersectRes is not None:
                        # print(intersectRes)

                        assert abs(intersectRes[0] - x) < 1e-6

                        node1 = [x, np.floor((intersectRes[1]))]
                        node2 = [x, np.ceil((intersectRes[1]))]

                        constraintId = [int(flatten2DIndex(*node1, s.gridSize)), int(flatten2DIndex(*node2, s.gridSize))]
                        weights = [node2[1] - intersectRes[1] , intersectRes[1] - node1[1]]

                        s.equalityConstraintIds.append(constraintId)
                        s.equalityConstraintWeights.append(weights)
                        s.equalityConstraintVal.append(inputContourLineData["ContourLineHeights"][iContour])
                for y in ys:
                    scanline = np.array([
                        [xMin - 1, y],
                        [xMax + 1, y],
                    ])

                    intersectRes = s.calculate_intersection(scanline, lineSeg)
                    if intersectRes is not None:
                        # print(intersectRes)

                        assert abs(intersectRes[1] - y) < 1e-6

                        node1 = [np.floor((intersectRes[0])), y]
                        node2 = [np.ceil((intersectRes[0])), y]

                        constraintId = [int(flatten2DIndex(*node1, s.gridSize)), int(flatten2DIndex(*node2, s.gridSize))]
                        weights = [node2[0] - intersectRes[0], intersectRes[0] - node1[0]]

                        s.equalityConstraintIds.append(constraintId)
                        s.equalityConstraintWeights.append(weights)
                        s.equalityConstraintVal.append(inputContourLineData["ContourLineHeights"][iContour])

        # plot
        if plot:
            pts = []
            for iConstraint in range(len(s.equalityConstraintIds)):
                constraintIds = s.equalityConstraintIds[iConstraint]
                weights = s.equalityConstraintWeights[iConstraint]
                constraintPoints = \
                    np.array([np.array(to2DIndex(i, s.gridSize) )*w for i, w in zip(constraintIds,weights )])

                p = np.sum(constraintPoints, axis=0)
                pts.append(p)

            pts = np.array(pts)
            s.scatter.set_offsets(np.c_[pts[:,0], pts[:,1]])
            s.fig.canvas.draw_idle()
            plt.pause(0.01)
            # plt.waitforbuttonpress()

    def saveConstraints(s, outFile, xRange, yRange):
        allPts = []
        for iConstraint in range(len(s.equalityConstraintIds)):
            constraintIds = s.equalityConstraintIds[iConstraint]
            weights = s.equalityConstraintWeights[iConstraint]
            constraintPoints = \
                np.array([np.array(to2DIndex(i, s.gridSize)) * w for i, w in zip(constraintIds, weights)])

            p = np.sum(constraintPoints, axis=0)
            x = xRange[0] + (xRange[1] - xRange[0]) * p[0] / s.gridSize[0]
            y = yRange[0] + (yRange[1] - yRange[0]) * p[1] / s.gridSize[1]

            allPts.append([x, y, s.equalityConstraintVal[iConstraint]])

        pc = pv.PolyData(allPts)
        pc.save(outFile, binary=False)
