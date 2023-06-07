import os

if __name__ == "__main__":

    for i in range(1, 100):
        t = i / 100
        print(t)

        input_dir = "./VortexSlice/"
        output_dir = "./VortexSlice/monoMesh_sub_8_10/"
        inputStartFile = "monoMesh_sub_008.vti"
        inputEndFile = "monoMesh_sub_010.vti"
        threshold = 0.02
        absoluteParam = 0  # 0 - relative persisththt7ht7hh7tht7ht7tht7tence threshold (1%) 1 - absolute persistence threshold (0.06)
        print(absoluteParam)
        treeType = "mt"
        scalarField = 'speed'

        submitCommand = "pvpython geodesic_other.py " + str(t) + " " + input_dir + \
                        " " + output_dir + " " + inputStartFile + " " + inputEndFile + " " + \
                        str(threshold) + " " + str(absoluteParam) + " " + treeType + " " + scalarField
        print(submitCommand)
        os.system(submitCommand)