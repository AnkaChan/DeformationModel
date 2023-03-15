import os

if __name__ == "__main__":
    for i in range(1, 100):
        t = i / 100
        print(t)

        submitCommand = "pvpython Data/geodesics/geodesic.py " + str(t)
        print(submitCommand)
        os.system(submitCommand)