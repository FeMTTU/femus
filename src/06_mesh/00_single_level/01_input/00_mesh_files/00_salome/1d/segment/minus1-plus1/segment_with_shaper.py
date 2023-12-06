#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/src/06_mesh/00_single_level/01_input/00_mesh_files/00_salome/01_1d/segment/minus1-plus1')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("n_x", 1)
####################################################
##        End of NoteBook variables section       ##
####################################################
###
### SHAPER component
###

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()

### Create Part
Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()
model.addParameter(Part_1_doc, "DeltaX", '2')
model.addParameter(Part_1_doc, "x_begin", '-1')
model.addParameter(Part_1_doc, "x_end", 'x_begin + DeltaX')

### Create Point
Point_2 = model.addPoint(Part_1_doc, "x_begin", "0", "0")

### Create Point
Point_3 = model.addPoint(Part_1_doc, "x_end", "0", "0")

### Create Edge
Edge_1 = model.addEdge(Part_1_doc, model.selection("VERTEX", "Point_1"), model.selection("VERTEX", "Point_2"))

### Create Export
Export_1 = model.exportToXAO(Part_1_doc, '/tmp/shaper_9vuagw9h.xao', model.selection("EDGE", "Edge_1_1"), 'XAO')

model.end()

###
### SHAPERSTUDY component
###

model.publishToShaperStudy()
import SHAPERSTUDY
Edge_1_1, = SHAPERSTUDY.shape(model.featureStringId(Edge_1))
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
(imported, Edge_1_1, [], [], []) = geompy.ImportXAO("/tmp/shaper_9vuagw9h.xao")
[Vertex_1,Vertex_2] = geompy.ExtractShapes(Edge_1_1, geompy.ShapeType["VERTEX"], True)
[Vertex_1, Vertex_2] = geompy.GetExistingSubObjects(Edge_1_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Edge_1_1, 'Edge_1_1' )
geompy.addToStudyInFather( Edge_1_1, Vertex_1, 'Vertex_1' )
geompy.addToStudyInFather( Edge_1_1, Vertex_2, 'Vertex_2' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Edge_1_1,'Mesh_1')
Regular_1D = Mesh_1.Segment()
Number_of_Segments_1 = Regular_1D.NumberOfSegments("n_x")
Group_1_0 = Mesh_1.GroupOnGeom(Vertex_1,'Vertex_1',SMESH.NODE)
Group_2_0 = Mesh_1.GroupOnGeom(Vertex_2,'Vertex_2',SMESH.NODE)
isDone = Mesh_1.Compute()
[ Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0)
Group_1_0.SetName( 'Group_1_0' )
Group_2_0.SetName( 'Group_2_0' )
[ Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
#try:
#  Mesh_1.ExportMED( r'/home/gbornia/software/femus/src/06_mesh/00_single_level/01_input/00_mesh_files/00_salome/01_1d/segment/minus1-plus1/segment_ndivisions_1.med', 0, 41, 1, Mesh_1, 0, [], '',-1, 1 )
#  pass
#except:
#  print('ExportPartToMED() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
