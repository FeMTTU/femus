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

Mesh_1 = smesh.Mesh(Face_2)
Regular_1D = Mesh_1.Segment()
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Number_of_Segments_1 = Regular_1D.NumberOfSegments(15)
isDone = Mesh_1.Compute()
Mesh_2 = smesh.Mesh(Face_3)
Regular_1D_1 = Mesh_2.Segment()
Number_of_Segments_2 = Regular_1D_1.NumberOfSegments(15)
Quadrangle_2D_1 = Mesh_2.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_2.Compute()
Mesh_3 = smesh.Mesh(Face_4)
Regular_1D_2 = Mesh_3.Segment()
Number_of_Segments_3 = Regular_1D_2.NumberOfSegments(15)
Quadrangle_2D_2 = Mesh_3.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_3.Compute()
Number_of_Segments_4 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_4.SetNumberOfSegments( 15 )
Mesh_4 = smesh.Mesh(Face_5)
Regular_1D_3 = Mesh_4.Segment()
Number_of_Segments_5 = Regular_1D_3.NumberOfSegments(15)
Quadrangle_2D_3 = Mesh_4.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_4.Compute()
Compound_Mesh_1 = smesh.Concatenate( [ Mesh_1.GetMesh(), Mesh_2.GetMesh(), Mesh_3.GetMesh(), Mesh_4.GetMesh() ], 1, 1, 1e-05, False )


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Number_of_Segments_5, 'Number of Segments_5')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Mesh_3.GetMesh(), 'Mesh_3')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')
smesh.SetName(Compound_Mesh_1.GetMesh(), 'Compound_Mesh_1')
smesh.SetName(Mesh_4.GetMesh(), 'Mesh_4')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
