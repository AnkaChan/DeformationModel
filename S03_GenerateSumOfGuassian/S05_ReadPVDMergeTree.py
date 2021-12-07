import pyvista as pv

if __name__ == '__main__':
    tree1File = r'X:\Code\CompTopology\Data\S04_ManuallyDesignMergeTree\MergeTree1.vtk'

    tree1 = pv.read(tree1File)
    print(tree1)
