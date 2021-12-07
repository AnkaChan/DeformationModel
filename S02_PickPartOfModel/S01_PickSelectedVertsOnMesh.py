from M02_ObjConverter import *

if __name__ == '__main__':
    inArmVertIdsFile = r'LeftArmVIds.json'
    inputFolder = r'F:\WorkingCopy2\2020_11_26_SMPLSHFit\Fit\Lada_Stand'
    outputFolder = r'F:\WorkingCopy2\2021_02_15_TestIsometricDeformation\Input\ArmJoint'
    exampleQuadMesh = 'SMPLWithSocks_Quad_xyzOnly_tri.obj'

    armVerts = json.load(open(inArmVertIdsFile))
    vertsToRemove = [i for i in range(6750) if i not in armVerts]

    removeVertsFromMeshFolder(inputFolder, outputFolder, vertsToRemove,
                              exampleQuadMesh, removeVerts=True, interval=[368, 408])