import os

if __name__ == "__main__":

    for i in range(1, 100):
        t = i / 100
        print(t)

        input_dir = "./HeatedFlowY/"
        output_dir = "./HeatedFlowY/data_1_3/"
        inputStartFile = "data_601.vti"
        inputEndFile = "data_603.vti"
        threshold = 0.1
        absoluteParam = 0  # 0 - relative persistence threshold (1%) 1 - absolute persistence threshold (0.06)
        print(absoluteParam)
        treeType = "st"
        scalarField = 'velocityMagnitude'

        submitCommand = "pvpython geodesic_other.py " + str(t) + " " + input_dir + \
                        " " + output_dir + " " + inputStartFile + " " + inputEndFile + " " + \
                        str(threshold) + " " + str(absoluteParam) + " " + treeType + " " + scalarField
        print(submitCommand)
        os.system(submitCommand)