#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/student/software/femus/applications/2021_fall/GayaniBalasuriya/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("r", 1)
notebook.set("x", 1)
notebook.set("y", 1)
####################################################
##        End of NoteBook variables section       ##
####################################################
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
Disk_2 = geompy.MakeDiskR("r", 1)
[geomObj_1] = geompy.SubShapeAll(Disk_2, geompy.ShapeType["EDGE"])
Face_1 = geompy.MakeFaceHW("x", "y", 1)
Vertex_1 = geompy.MakeVertex(0, 0, 0)
Face_1_vertex_9 = geompy.GetSubShape(Face_1, [9])
Face_1_vertex_5 = geompy.GetSubShape(Face_1, [5])
Face_1_vertex_4 = geompy.GetSubShape(Face_1, [4])
Face_1_vertex_7 = geompy.GetSubShape(Face_1, [7])
[Face_1_vertex_9, Face_1_vertex_5, Face_1_vertex_4, Face_1_vertex_7] = geompy.GetExistingSubObjects(Face_1, False)
[geomObj_1] = geompy.GetExistingSubObjects(Disk_2, False)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Face_1_vertex_9)
Scale_1 = geompy.MakeScaleTransform(Line_1, None, 2)
Scale_1_vertex_2 = geompy.GetSubShape(Scale_1, [2])
Line_2 = geompy.MakeLineTwoPnt(Scale_1_vertex_2, Face_1_vertex_5)
Scale_2 = geompy.MakeScaleTransform(Line_2, None, 2)
Scale_2_vertex_2 = geompy.GetSubShape(Scale_2, [2])
Line_3 = geompy.MakeLineTwoPnt(Scale_2_vertex_2, Face_1_vertex_4)
Scale_3 = geompy.MakeScaleTransform(Line_3, None, 1.5)
Scale_3_vertex_2 = geompy.GetSubShape(Scale_3, [2])
Line_4 = geompy.MakeLineTwoPnt(Scale_3_vertex_2, Face_1_vertex_7)
Scale_4 = geompy.MakeScaleTransform(Line_4, None, 1.5)
Disk_2_edge_3 = geompy.GetSubShape(Disk_2, [3])
Vertex_2 = geompy.MakeVertexOnLinesIntersection(Scale_3, Disk_2_edge_3)
Vertex_3 = geompy.MakeVertexOnLinesIntersection(Scale_2, Disk_2_edge_3)
Vertex_4 = geompy.MakeVertexOnLinesIntersection(Scale_1, Disk_2_edge_3)
Vertex_5 = geompy.MakeVertexOnLinesIntersection(Scale_4, Disk_2_edge_3)
Line_5 = geompy.MakeLineTwoPnt(Face_1_vertex_5, Vertex_3)
Line_6 = geompy.MakeLineTwoPnt(Face_1_vertex_4, Vertex_2)
Line_7 = geompy.MakeLineTwoPnt(Face_1_vertex_7, Vertex_5)
Line_8 = geompy.MakeLineTwoPnt(Face_1_vertex_9, Vertex_4)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Disk_2, 'Disk_2' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudyInFather( Face_1, Face_1_vertex_9, 'Face_1:vertex_9' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Scale_1, 'Scale_1' )
geompy.addToStudyInFather( Scale_1, Scale_1_vertex_2, 'Scale_1:vertex_2' )
geompy.addToStudyInFather( Face_1, Face_1_vertex_5, 'Face_1:vertex_5' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Scale_2, 'Scale_2' )
geompy.addToStudyInFather( Scale_2, Scale_2_vertex_2, 'Scale_2:vertex_2' )
geompy.addToStudyInFather( Face_1, Face_1_vertex_4, 'Face_1:vertex_4' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Scale_3, 'Scale_3' )
geompy.addToStudyInFather( Scale_3, Scale_3_vertex_2, 'Scale_3:vertex_2' )
geompy.addToStudyInFather( Face_1, Face_1_vertex_7, 'Face_1:vertex_7' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudy( Scale_4, 'Scale_4' )
geompy.addToStudyInFather( Disk_2, Disk_2_edge_3, 'Disk_2:edge_3' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudy( Line_6, 'Line_6' )
geompy.addToStudy( Line_7, 'Line_7' )
geompy.addToStudy( Line_8, 'Line_8' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Disk_2)
Regular_1D = Mesh_1.Segment()
#Group_1 = Mesh_1.GroupOnGeom(__NOT__Published__Object__,'',SMESH.EDGE)
Mesh_2 = smesh.Mesh(Face_1)
Quadrangle_2D = Mesh_2.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Face_1_vertex_9_1 = Mesh_2.GroupOnGeom(Face_1_vertex_9,'Face_1:vertex_9',SMESH.NODE)
Face_1_vertex_5_1 = Mesh_2.GroupOnGeom(Face_1_vertex_5,'Face_1:vertex_5',SMESH.NODE)
Face_1_vertex_4_1 = Mesh_2.GroupOnGeom(Face_1_vertex_4,'Face_1:vertex_4',SMESH.NODE)
Face_1_vertex_7_1 = Mesh_2.GroupOnGeom(Face_1_vertex_7,'Face_1:vertex_7',SMESH.NODE)


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Face_1_vertex_7_1, 'Face_1:vertex_7')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Face_1_vertex_9_1, 'Face_1:vertex_9')
smesh.SetName(Face_1_vertex_4_1, 'Face_1:vertex_4')
smesh.SetName(Face_1_vertex_5_1, 'Face_1:vertex_5')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
