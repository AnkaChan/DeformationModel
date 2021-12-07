from M01_TopologicalExtraction import *
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import cv2

if __name__ == '__main__':
    inputScalarField = r'X:\Data\01_Topologgy\DatasetsVTKFormat\DatasetsVTKFormat\HeatedFlow\TopologicalData\data_600\data_600.vti'
    inputSegmentedFieldSource = r'X:\Data\01_Topologgy\DatasetsVTKFormat\DatasetsVTKFormat\HeatedFlow\TopologicalData\data_600\Seg.vti'
    inputSegmentedFieldTarget = r'X:\Data\01_Topologgy\DatasetsVTKFormat\DatasetsVTKFormat\HeatedFlow\TopologicalData\data_630\Seg.vti'
    inputMergeTreeNodesField = r'X:\Code\CompTopology\Data\S04_GenerateNewScalarField.py\Topology\Node.vtk'
    outFolder = r'F:\WorkingCopy2\2021_07_13_ContourLineDeformation\HeatMap'

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(inputSegmentedFieldSource)

    reader.Update()
    image = reader.GetOutput()

    # # print(image)
    #
    # dims = image.GetDimensions()
    # # dims = (width, height) = (cols, rows)
    # dimsRC = (dims[1], dims[0])
    #
    # scalar = vtk_to_numpy(image.GetPointData().GetArray("SegmentationId"))
    #
    # for iRow in range(dimsRC[0]):
    #     scalar[flatten2DIndex(iRow, int(dimsRC[1]/2), dimsRC), ] = 0
    # segImg = scalar.reshape(dimsRC)
    # # ‘C’ means to read / write the elements using C-like index order,
    # # with the last axis index changing fastest, back to the first axis index changing slowest.
    # # such scalar is row major
    #
    # segImg = 255*segImg/np.max(segImg)
    #
    # cv2.imwrite("seg.png", np.squeeze(segImg.astype(np.uint8)))

    # how to process vti file?
    # how to read its attibutes?

    segs, dimsRC = getVTIData(image, "SegmentationId", reshape2D=True)
    initalLSF = np.zeros(dimsRC, dtype=np.uint8)
    initalLSF[np.logical_or(segs==43, segs==44)] = 255

    cv2.imshow('LSF', initalLSF)
    cv2.waitKey()

    segs, dimsRC =  readVTIData(inputSegmentedFieldTarget, "SegmentationId", reshape2D=True)

    targetLSF = np.zeros(dimsRC, dtype=np.uint8)
    targetLSF[np.logical_or(segs==43, segs==44)] = 255

    cv2.imshow('LSF', targetLSF)
    cv2.waitKey()

    cv2.imwrite(join(outFolder, 'Source.png'), initalLSF)
    cv2.imwrite(join(outFolder, 'Target.png'), targetLSF)

    # target =


    # contourEdges, contourLineConstraintWeight, contourLineConstraintHeight = findContourLineConstraints(segs, dimsRC, contourLineHeight)


