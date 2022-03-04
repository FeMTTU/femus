#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/student/software/femus/applications/2021_fall/hdongamm/input')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Divided_Disk_1 = geompy.MakeDividedDisk(1, 1, GEOM.SQUARE)
geompy.Rotate(Divided_Disk_1, OZ, 45*math.pi/180.0)
Face_1 = geompy.MakeFaceHW(2, 2, 1)
Translation_1 = geompy.MakeTranslation(Face_1, 0, -1, 0)
Cut_1 = geompy.MakeCutList(Divided_Disk_1, [Translation_1], True)
[Face_2,Face_3,Face_4,Face_5] = geompy.ExtractShapes(Cut_1, geompy.ShapeType["FACE"], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Divided_Disk_1, 'Divided Disk_1' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudyInFather( Cut_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Cut_1, Face_3, 'Face_3' )
geompy.addToStudyInFather( Cut_1, Face_4, 'Face_4' )
geompy.addToStudyInFather( Cut_1, Face_5, 'Face_5' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Regular_1D_4 = smesh.CreateHypothesis( "Regular_1D" )
Quadrangle_2D_4 = smesh.CreateHypothesis( "Quadrangle_2D" )
Mesh_2 = smesh.Mesh(Face_3)
#status = Mesh_2.AddHypothesis(Number_of_Segments_1) ### Number_of_Segments_1 has not been yet created
Regular_1D_1 = Mesh_2.Segment()
Quadrangle_2D_1 = Mesh_2.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_2.Compute()
Mesh_3 = smesh.Mesh(Face_4)
#status = Mesh_3.AddHypothesis(Number_of_Segments_1) ### Number_of_Segments_1 has not been yet created
Regular_1D_2 = Mesh_3.Segment()
Quadrangle_2D_2 = Mesh_3.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_3.Compute()
Mesh_4 = smesh.Mesh(Face_5)
#status = Mesh_4.AddHypothesis(Number_of_Segments_1) ### Number_of_Segments_1 has not been yet created
Regular_1D_3 = Mesh_4.Segment()
Quadrangle_2D_3 = Mesh_4.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_4.Compute()
#Mesh_1 = smesh.Concatenate( [ Mesh_2.GetMesh() ], 1, 1, 1e-05, True ,meshToAppendTo=Mesh_1.GetMesh()) ### Mesh_1 has not been yet created
#[ GrMesh_2_Nodes, GrMesh_2_Edges, GrMesh_2_Faces ] = Mesh_1.GetGroups() ### not created Object
#Compound_Mesh_1 = smesh.Concatenate( [ Mesh_1.GetMesh(), Mesh_2.GetMesh(), Mesh_3.GetMesh(), Mesh_4.GetMesh() ], 1, 1, 1e-05, False ) ### Mesh_1 has not been yet created
#[ GrMesh_2_Nodes_1, GrMesh_2_Edges_1, GrMesh_2_Faces_1 ] = Compound_Mesh_1.GetGroups() ### not created Object


## Set names of Mesh objects
smesh.SetName(Mesh_4.GetMesh(), 'Mesh_4')
smesh.SetName(Mesh_3.GetMesh(), 'Mesh_3')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
