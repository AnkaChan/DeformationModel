import MeshFrame as MF
import numpy as np
mesh = MF.Mesh()
mesh.load('example.obj')

vNormalHanle = mesh.addVProp("Normal")

for vId in range(len(mesh.vertices())):
    totalArea = 0
    vNormal = np.zeros((3,), dtype=np.float64)
    for fId in mesh.VFIter(vId):
        fNormal = mesh.faces(fId).normal()
        fArea = mesh.faces(fId).area()
        vNormal = vNormal + fArea * fArea
        totalArea += fArea

    mesh.setVProp(vId, vNormalHanle, vNormal / totalArea)












