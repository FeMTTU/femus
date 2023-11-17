#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/src/06_mesh/00_single_level/01_input/00_mesh_files/00_salome/01_1d/segment/0-1')




###
### SHAPER component
###

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()

### Create Part
Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()

### Create Point
Point_2 = model.addPoint(Part_1_doc, 0, 0, 0)

### Create Point
Point_3 = model.addPoint(Part_1_doc, 1, 0, 0)

### Create Sketch
Sketch_1 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_1 = Sketch_1.addLine(0, 0, 1, 0)

### Create SketchProjection
SketchProjection_1 = Sketch_1.addProjection(model.selection("VERTEX", "Point_1"), False)
SketchPoint_1 = SketchProjection_1.createdFeature()
Sketch_1.setCoincident(SketchLine_1.startPoint(), SketchPoint_1.result())

### Create SketchProjection
SketchProjection_2 = Sketch_1.addProjection(model.selection("VERTEX", "Point_2"), False)
SketchPoint_2 = SketchProjection_2.createdFeature()
Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchPoint_2.result())
model.do()


### Create Vertex
Vertex_1 = model.addVertex(Part_1_doc, [model.selection("VERTEX", "all-in-Point_1")], False)

### Create Vertex
Vertex_2 = model.addVertex(Part_1_doc, [model.selection("VERTEX", "Point_2")], False)

### Create Edge
Edge_1 = model.addEdge(Part_1_doc, [model.selection("EDGE", "Sketch_1/SketchLine_1")])

### Create Export
Export_3 = model.exportToXAO(Part_1_doc, '/tmp/shaper_z2369cmc.xao', model.selection("VERTEX", "Vertex_1_1"), 'XAO')

### Create Export
Export_4 = model.exportToXAO(Part_1_doc, '/tmp/shaper_dfclz5n1.xao', model.selection("VERTEX", "Vertex_2_1"), 'XAO')

### Create Export
Export_5 = model.exportToXAO(Part_1_doc, '/tmp/shaper_9u78ow0x.xao', model.selection("EDGE", "Edge_1_1"), 'XAO')

model.end()

###
### SHAPERSTUDY component
###

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
(imported, Vertex_1_1, [], [], []) = geompy.ImportXAO("/tmp/shaper_z2369cmc.xao")
(imported, Vertex_2_1, [], [], []) = geompy.ImportXAO("/tmp/shaper_dfclz5n1.xao")
(imported, Edge_1_1, [], [], []) = geompy.ImportXAO("/tmp/shaper_9u78ow0x.xao")
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Vertex_1_1, 'Vertex_1_1' )
geompy.addToStudy( Vertex_2_1, 'Vertex_2_1' )
geompy.addToStudy( Edge_1_1, 'Edge_1_1' )

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
Number_of_Segments_1 = Regular_1D.NumberOfSegments(1)  # "n_x"
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0)
smesh.SetName(Mesh_1, 'Mesh_1')

smeshObj_1 = Mesh_1.CreateEmptyGroup( SMESH.NODE, 'Group_1' )
nbAdd = smeshObj_1.Add( [ 1, 2 ] )
smeshObj_1.SetName( 'Group_1_0' )
[ smeshObj_1 ] = Mesh_1.GetGroups()
Mesh_1.RemoveGroup( smeshObj_1 )
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.NODE, 'Group_1_0' )
nbAdd = Group_1_0.Add( [ 1 ] )
Group_2_1 = Mesh_1.CreateEmptyGroup( SMESH.NODE, 'Group_2_1' )
nbAdd = Group_2_1.Add( [ 2 ] )

## some objects were removed
aStudyBuilder = salome.myStudy.NewBuilder()
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_1))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)

## Set names of Mesh objects
smesh.SetName(Group_2_1, 'Group_2_1')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
