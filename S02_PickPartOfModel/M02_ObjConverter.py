import glob
import os, tqdm
from pathlib import Path
import pyvista as pv
from os.path import join
import json
from SkelFit.Data import retrieveFaceStructure, flattenFaces
import numpy as np

# vt_path = r'..\Data\BuildSmplsh\SMPLWithSocks_tri.obj'
# vt_path = r'..\Data\BuildSmplsh\SMPLWithSocks_Quad.obj'
#
# vts = []
# fs = []
# vns = []
# vs = []
# with open(vt_path, 'r') as f:
#     lines = f.readlines()
#     for line in lines:
#         l = line.split(' ')
#         if l[0] == 'vt':
#             u = l[1]
#             v = l[2].split('\n')[0]
#             vts.append([u, v])
#         elif l[0] == 'vn':
#             vns.append([0, 0, 0])
#         elif l[0] == 'v':
#             vs.append([0, 0, 0])
#         elif l[0] == 'f':
#             fs_curr = []
#             for i in range(len(l) - 1):
#                 fi = l[i+1].split('/')
#                 fi = '{}/{}/{}'.format(fi[0], fi[1], fi[2].split('\n')[0])
#                 fs_curr.append(fi)
#             fs.append(fs_curr)
#     f.close()
# fs = json.load(open('SMPLSHQuadFaces.json'))
# print(len(vns))
# print(len(vs))
# print(len(vts))
# print(len(fs))


def converObjsInFolder(obj_dir, out_dir, ext='obj', convertToMM=False, addA=False, withMtl = False,  textureFile='', facesOnSuit=None):
    os.makedirs(out_dir, exist_ok=True)

    in_paths = glob.glob(obj_dir + '/*.' + ext)
    for in_path in in_paths:
        obj_name = Path(in_path).stem
        if addA:
            out_path = out_dir + '/A{}.obj'.format(obj_name)
        else:
            out_path = out_dir + '/{}.obj'.format(obj_name)

        convertObjFile(in_path, out_path, convertToMM, withMtl=withMtl, textureFile=textureFile, facesToPreserve=facesOnSuit)


def removeVertsFromMeshFolder(inFolder,  outFolder, vertIdsToRemove, exampleMesh, inExtName='ply', outExtName='ply',  removeVerts=True, interval=None):
    os.makedirs(outFolder, exist_ok=True)

    vertIdsToRemove = set(vertIdsToRemove)
    exampleMesh = pv.PolyData(exampleMesh)
    faces = retrieveFaceStructure(exampleMesh)
    vertIds = list(range(exampleMesh.points.shape[0]))

    if removeVerts:
        for iV in vertIdsToRemove:
            vertIds.remove(iV)

        vertToOldId = {vId: iV for iV, vId in enumerate(vertIds)}
        faces = retrieveFaceStructure(exampleMesh)
        newFaces = [[vertToOldId[iFV] for iFV in faces[iF]] for iF in range(len(faces)) if
                    np.all([iV not in vertIdsToRemove for iV in faces[iF]])]
    else:
        newFaces = [faces[iF] for iF in range(len(faces)) if np.all([iV not in vertIdsToRemove for iV in faces[iF]])]

    inFiles = glob.glob(join(inFolder, '*.'+inExtName))
    if interval is not None:
        inFiles = inFiles[interval[0]: interval[1]]

    for inFile in tqdm.tqdm(inFiles, desc='Removing Faces'):
        fileName = Path(inFile).stem
        outFile = join(outFolder, fileName+'.'+outExtName)

        removeVertsFromMesh(inFile, outFile, vertIds, exampleMesh, newFaces=newFaces)

def removeVertsFromMesh(inMesh, outMesh, vertIdsToKeep, exampleMesh, newFaces=None, removeVerts=True):
    inMesh = pv.PolyData(inMesh)

    # newVerts = [inMesh.points[iV] for iV in range(inMesh.points.shape[0]) if iV not in vertIdsToRemove]
    if removeVerts:
        newVerts = inMesh.points[vertIdsToKeep, :].tolist()
    else:
        newVerts = [inMesh.points[iV] for iV in range(inMesh.points.shape[0]) ]

    if newFaces is not None:
        newFacesFlattern = flattenFaces(newFaces)
    else:
        exampleMesh = pv.PolyData(exampleMesh)
        faces = retrieveFaceStructure(exampleMesh)
        vertToOldId = {vId: iV for iV, vId in enumerate(vertIdsToKeep)}
        newFaces = [[vertToOldId[iFV] for iFV in faces[iF]] for iF in range(len(faces)) if np.all([iV in vertIdsToKeep for iV in faces[iF]])]
        newFacesFlattern = flattenFaces(newFaces)

    newMesh = pv.PolyData(np.array(newVerts), np.array(newFacesFlattern))

    newMesh.save(outMesh)

def convertObjFile(inFile, outFile, convertToMM=False, withMtl=False, textureFile=None, facesToPreserve=None):
    # obj_name = in_path.split('\\')[-1]

    extName = Path(inFile).suffix
    # read current
    vs = []

    if extName.lower() == '.obj':
        with open(inFile, 'r') as f:
            lines = f.readlines()

            for line in lines:
                l = line.split(' ')
                if l[0] == 'v':
                    vs.append([l[1], l[2], l[3].split('\n')[0]])
                elif l[0] == 'vt' or l[0] == 'vn':
                    continue
            f.close()
    else:
        mesh = pv.PolyData(inFile)
        vs = mesh.points.tolist()

    # write new
    with open(outFile, 'w+') as f:
        fp = Path(outFile)
        outMtlFile = join(str(fp.parent), fp.stem + '.mtl')
        if withMtl:
            f.write('mtllib ./' + fp.stem + '.mtl\n')
            with open(outMtlFile, 'w') as fMtl:
                mtlStr = '''newmtl material_0
Ka 0.200000 0.200000 0.200000
Kd 1.000000 1.000000 1.000000
Ks 1.000000 1.000000 1.000000
Tr 1.000000
illum 2
Ns 0.000000
map_Kd '''
                mtlStr += textureFile
                fMtl.write(mtlStr)


        for i, v in enumerate(vs):
            vn = vns[i]

            f.write('vn {} {} {}\n'.format(vn[0], vn[1], vn[2]))
            if convertToMM:
                v[0] = 1000*v[0]
                v[1] = 1000*v[1]
                v[2] = 1000*v[2]
            f.write('v {} {} {}\n'.format(v[0], v[1], v[2]))
        for vt in vts:
            f.write('vt {} {}\n'.format(vt[0], vt[1]))

        if withMtl:
            f.write('usemtl material_0\n')
        for iF in range(len(fs)):
            if facesToPreserve is not None and iF not in facesToPreserve:
                continue
            f.write('f')
            for fi in fs[iF]:
                f.write(' {}'.format(fi))
            f.write('\n')
        f.close()
    print(outFile)
    print(len(vs))

def objFilesToPly(inFolder, outFolder):
    os.makedirs(outFolder, exist_ok=True)
    objFiles = glob.glob(join(inFolder, '*.obj'))

    for objFile in objFiles:
        deformedMesh = pv.PolyData(objFile)
        outFile = join(outFolder, os.path.basename(objFile) + '.ply')
        deformedMesh.save(outFile, binary=False)

if __name__ == '__main__':
    # obj_dir = r'F:\WorkingCopy2\2020_07_15_NewInitialFitting\IniitalTexture\Meshes'
    # out_dir = r'E:\WorkingCopy\2020_06_30_AC_ConsequtiveTexturedFitting2\FinalObj\WithTextureCoord'
    # obj_dir = r'F:\WorkingCopy2\2020_07_28_TexturedFitting_Lada\Final\Mesh'
    # obj_dir = r'F:\WorkingCopy2\2020_08_27_KateyBodyModel\InitialSilhouetteFitting_NoGlassese\Final\18411'
    # obj_dir = r'C:\Code\MyRepo\03_capture\Mocap-CVPR-Paper-Figures\09_PipelineALL\Data\ObjWithTexture'
    # obj_dir = r'C:\Code\MyRepo\03_capture\BodyTracking\Data\KateyBodyModel\BodyMesh\Initial'
    # obj_dir = r'F:\WorkingCopy2\2020_08_26_TexturedFitting_LadaGround\FitOnlyBody\Vis\ObjWithUV'
    # obj_dir = r'C:\Code\MyRepo\03_capture\Mocap-CVPR-Paper-Figures\12_TeaserImage\Mesh'

    # obj_dir = r'F:\WorkingCopy2\2020_07_15_NewInitialFitting\CompleteTexture\Meshes'
    # texture = r'texturemap_learned_LapW0.2_MaskTrue_L1.png'
    #
    # objFilesToPly(r'E:\Dropbox\mcproj\2021_01_05_MoreAnimationSeqs\Katey', join(r'E:\Dropbox\mcproj\2021_01_05_MoreAnimationSeqs\Katey', 'ply'))
    # objFilesToPly(r'E:\Dropbox\mcproj\2021_01_05_MoreAnimationSeqs\LadaGround', join(r'E:\Dropbox\mcproj\2021_01_05_MoreAnimationSeqs\LadaGround', 'ply'))
    # objFilesToPly(r'E:\Dropbox\mcproj\2021_01_05_MoreAnimationSeqs\LadaStand', join(r'E:\Dropbox\mcproj\2021_01_05_MoreAnimationSeqs\LadaStand', 'ply'))

    # facesFile = 'FacesOnlySuit.json'
    # facesOnSuit = set(json.load(open(facesFile)))
    # facesOnSuit = None
    #
    # ext = 'obj'
    # # withMtl = True
    # withMtl = False
    # rename = False
    # plyFiles = glob.glob(join(obj_dir, '*.' + ext))
    # if rename:
    #     for plyF in plyFiles:
    #         fName = os.path.basename(plyF)
    #         if fName[0] != 'A':
    #             newFileName = join(obj_dir, 'A'+fName)
    #             os.rename(plyF, newFileName)
    #             # print(plyF, newFileName)
    #
    # out_dir = os.path.join(obj_dir, 'WithTextureCoord')
    # converObjsInFolder(obj_dir, out_dir, ext=ext, addA=False, withMtl=withMtl, textureFile=texture, facesOnSuit=facesOnSuit)
    #
    # objFilesToPly(out_dir, join(obj_dir, 'PlyWithTextureCoord'))

    # inFile = r'F:\WorkingCopy2\2020_07_15_NewInitialFitting\InitialSilhouetteFitting\3052\Final\InterpolatedWithSparse.ply'
    # outFile = r'F:\WorkingCopy2\2020_07_15_NewInitialFitting\InitialSilhouetteFitting\3052\Final\FinalMesh.obj'
    # convertObjFile(inFile, outFile)

    from SkelFit.Data import getIsolatedVerts
    # inFolder = r'F:\WorkingCopy2\2021_01_04_NewModelFitting\Output\Lada_Stand\SLap_SBiLap_True_TLap_1_JTW_0.5_JBiLap_0_Step1_Overlap0\Interpolated\Test'
    # inFolder = r'F:\WorkingCopy2\2021_01_04_NewModelFitting\Output\Katey_Interpolation\TLap_100_JR_600_JTW_200_Step150_Overlap50_cleaned\RestPoseTarget\InterpolationDisplacement'
    # outFolder = join(inFolder, 'RestPoseInterpolated')

    # inFolder = r'F:\WorkingCopy2\2021_01_14_AnimatinoSeqs\LongSequences\LadaGround'
    # outFolder = r'F:\WorkingCopy2\2021_01_21_DataToSubmit\FinalResult\Male_1'

    inFolder = r'F:\WorkingCopy2\2021_01_14_AnimatinoSeqs\LongSequences\LadaStand'
    outFolder = r'F:\WorkingCopy2\2021_01_21_DataToSubmit\FinalResult\Male_2'

    headVIdsFile = r'..\Data\2020_12_27_betterCoarseMesh\Mesh1487\HeadVIdsWithNeck.Json'
    handVIdsFile = r'..\Data\2020_12_27_betterCoarseMesh\Mesh1487\HandVIds.json'
    exampleQuadMesh = r'..\Data\2020_12_27_betterCoarseMesh\Mesh1487\Complete_withHeadHand_XYZOnly.obj'

    headVIds = json.load(open(headVIdsFile))
    handVIds = json.load(open(handVIdsFile))
    isolatedPoints = getIsolatedVerts(pv.PolyData(exampleQuadMesh)).tolist()
    vertsToRemove = set(headVIds + handVIds + isolatedPoints)

    removeVertsFromMeshFolder(join(inFolder), outFolder, vertsToRemove, exampleQuadMesh, removeVerts=True)