from M01_TopologicalExtraction import *

import pyvista as pv
from os.path import join
import vtk
import numpy as np
from qpsolvers import solve_qp, available_solvers
from scipy import sparse

def interpolateScalarField(gridSize, P, equalityConstraintIds, equalityConstraintWeights, equalityConstraintVal,
                           inequalityConstraintIds, inequalityConstraintWeights, inEqualityConstraintVal, fixBoundary=False, boundaryValue=0,
                           dType=np.float64, verbose=True, solver='osqp'):
    '''
    convert equality and inequality constraints into qradratic programming form and solve it.
    gridSize: dimensions of the scalar field
    P: the image laplacian matrix of the input scalar field, s.t., X^T A X is the image laplacian of the input scalar field where
        X is all the values of the input scalar field. Must be a sparse.csc_matrix.
    equalityConstraintIds: a list contains N lists, N is the number of equality constraints. Each sub-list in equalityConstraintIds has
        different lengths, and it contains the flatten index of the pixels, that are under this equality constraint.
    equalityConstraintWeights: a list contains N lists, N is the number of equality constraints. Each sub-list in equalityConstraintWeights has
        different lengths, and its length is equal to the length of corresponding sub-list in equalityConstraintIds. Each sub-list contains
        floating point numbers and they sum up to 1 if it is a contour line or extremity constraint.
    equalityConstraintVal: a list contains N floating point numbers, N is the number of equality constraints.
        Say a sub-list in equalityConstraintIds is [id1, id2, id3, ...], its corresponding sub-list in equalityConstraintWeights
        is [w1, w2, w3, ...], its corresponding value in equalityConstraintVal isv; then the containt will be:
        x_id1 * w1 + x_id2 * w2 + x_id3 * w3 + ... = v

    inequalityConstraintIds: similar to equalityConstraintIds, but for inequality constraints.
    inequalityConstraintWeights: similar to equalityConstraintWeights, but for inequality constraints and don't need to sum up to 1.
    inEqualityConstraintVal: equalityConstraintVal to equalityConstraintWeights, but for inequality constraints.
        Say a sub-list in inequalityConstraintIds is [id1, id2, id3, ...], its corresponding sub-list in inequalityConstraintWeights
        is [w1, w2, w3, ...], its corresponding value in inEqualityConstraintVal v; then the containt will be:
        x_id1 * w1 + x_id2 * w2 + x_id3 * w3 + ... <= v
    fixBoundary: whether to fix the (image) boundary values
    boundaryValue: the value that (image) boundary fix to
    dType: data type the solver will be using
    verbose: where output detailed solving information
    solver: type of solvers, please refer qpsolvers.solve_qp and qpsolvers.available_solvers for available options
    '''

    numberOfViables = P.shape[0]

    if fixBoundary:
        # get boundary vertices
        boundaryVerts = []
        for i in range(gridSize[0]):
            boundaryVerts.append([flatten2DIndex(i, 0, gridSize)])
            boundaryVerts.append([flatten2DIndex(i, gridSize[1] - 1, gridSize)])

        for j in range(1, gridSize[1] - 1):
            boundaryVerts.append([flatten2DIndex(0, j, gridSize)])
            boundaryVerts.append([flatten2DIndex(gridSize[0] - 1, j, gridSize)])

        print("Number of boundary constraints", len(boundaryVerts))

        # equality constraints
        equalityConstraintIds = boundaryVerts + equalityConstraintIds
        equalityConstraintWeights = [[1.0] for i in range(len(boundaryVerts))] + equalityConstraintWeights
        equalityConstraintVal = [boundaryValue for i in range(len(boundaryVerts))] + equalityConstraintVal

    A = np.zeros((len(equalityConstraintIds), numberOfViables)).astype(dType)
    for i, constraintVId in enumerate(equalityConstraintIds):
        assert len(equalityConstraintWeights[i]) == len(equalityConstraintIds[i])
        for iPixel in range(len(constraintVId)):
            A[i, int(constraintVId[iPixel])] = equalityConstraintWeights[i][iPixel]
    b = np.array(equalityConstraintVal).astype(dType)

    # inequality constraint
    G = np.zeros((len(inequalityConstraintIds), numberOfViables)).astype(dType)
    for i, constraintVId in enumerate(inequalityConstraintIds):
        assert len(inequalityConstraintWeights[i]) == len(inequalityConstraintIds[i])
        for iPixel in range(len(constraintVId)):
            G[i, int(constraintVId[iPixel])] = inequalityConstraintWeights[i][iPixel]

    h = np.array(inEqualityConstraintVal).astype(dType)

    q = np.zeros((numberOfViables,)).astype(dType)

    G = sparse.csc_matrix(G)
    A = sparse.csc_matrix(A)

    Z = solve_qp(P, q, G, h, A, b, verbose=verbose, solver=solver)

    Z = Z.reshape(gridSize)

    return Z

