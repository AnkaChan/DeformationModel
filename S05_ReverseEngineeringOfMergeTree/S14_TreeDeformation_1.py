import glob
import csv
import pandas as pd

from M01_TopologicalExtraction import *

if __name__ == '__main__':
    inDeformDataFolder = r'F:\Code\02_Graphics\DeformationModel\S05_ReverseEngineeringOfMergeTree\Data\geodesics\JulienExample\animation'

    nameToTree2s = "interpolatedTotree2"
    deformationDataFileToTree2s = glob.glob(join(inDeformDataFolder, nameToTree2s + "*.csv"))

    for stateFile in deformationDataFileToTree2s:
        file = open(stateFile)
        # rows = []
        # csvreader = csv.reader(file)
        # header = next(csvreader)
        # for row in csvreader:
        #     rows.append(row)
        # print(header)
        # print(rows)





