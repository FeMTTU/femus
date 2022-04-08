#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/student')

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
Disk_1 = geompy.MakeDiskR(1, 1)
Face_1 = geompy.MakeFaceHW(2, 2, 1)
Translation_1 = geompy.MakeTranslation(Face_1, 0, -1, 0)
Cut_1 = geompy.MakeCutList(Disk_1, [Translation_1], True)
[Edge_1,Edge_2] = geompy.ExtractShapes(Cut_1, geompy.ShapeType["EDGE"], True)
[Edge_1, Edge_2] = geompy.GetExistingSubObjects(Cut_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Disk_1, 'Disk_1' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudyInFather( Cut_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Cut_1, Edge_2, 'Edge_2' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Cut_1)
Regular_1D = Mesh_1.Segment()
Number_of_Segments_1 = Regular_1D.NumberOfSegments(15)
NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_2D)
Edge_1_1 = Mesh_1.GroupOnGeom(Edge_1,'Edge_1',SMESH.EDGE)
Edge_2_1 = Mesh_1.GroupOnGeom(Edge_2,'Edge_2',SMESH.EDGE)
isDone = Mesh_1.Compute()
[ Edge_1_1, Edge_2_1 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0)
[ Edge_1_1, Edge_2_1 ] = Mesh_1.GetGroups()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Edge_1_1, 'Edge_1')
smesh.SetName(Edge_2_1, 'Edge_2')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
