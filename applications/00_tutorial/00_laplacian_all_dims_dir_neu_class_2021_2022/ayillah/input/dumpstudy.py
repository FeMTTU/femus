#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/student/software/femus/applications/2021_fall/ayillah/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("l_x", 3)
notebook.set("l_y", 1)
notebook.set("l_x_half", "0.5*l_x")
notebook.set("l_y_half", "0.5*l_y")
notebook.set("l_x", 3)
notebook.set("l_y", 1)
notebook.set("l_x_half", "0.5*l_x")
notebook.set("l_y_half", "0.5*l_y")
notebook.set("n_x", 4)
notebook.set("n_y", 3)
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
Face_1 = geompy.MakeFaceHW("l_x", "l_y", 1)
[geomObj_1,geomObj_2,geomObj_3,geomObj_4] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)
[geomObj_5,geomObj_6,geomObj_7,geomObj_8] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)
Translation_1 = geompy.MakeTranslation(Face_1, "l_x_half", "l_y_half", 0)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Translation_1, geompy.ShapeType["EDGE"], True)
[Edge_1, Edge_2, Edge_3, Edge_4] = geompy.GetExistingSubObjects(Translation_1, False)
Vertex_1 = geompy.MakeVertex(0, 5, 0)
Vertex_2 = geompy.MakeVertex(7, 0, 0)
Vertex_3 = geompy.MakeVertex(7, 5, 0)
Vertex_4 = geompy.MakeVertex(0, 0, 0)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudyInFather( Translation_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Translation_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Translation_1, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Translation_1, Edge_4, 'Edge_4' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Translation_1)
Edge_1_1 = Mesh_1.GroupOnGeom(Edge_1,'Edge_1',SMESH.EDGE)
Edge_2_1 = Mesh_1.GroupOnGeom(Edge_2,'Edge_2',SMESH.EDGE)
Edge_3_1 = Mesh_1.GroupOnGeom(Edge_3,'Edge_3',SMESH.EDGE)
Edge_4_1 = Mesh_1.GroupOnGeom(Edge_4,'Edge_4',SMESH.EDGE)
Regular_1D = Mesh_1.Segment(geom=Edge_2)
Number_of_Segments_1 = Regular_1D.NumberOfSegments("n_x")
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Regular_1D_1 = Mesh_1.Segment(geom=Edge_1)
Number_of_Segments_2 = Regular_1D_1.NumberOfSegments("n_y")
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_1)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_2 = Mesh_1.Segment()
isDone = Mesh_1.Compute()
[ Edge_1_1, Edge_2_1, Edge_3_1, Edge_4_1 ] = Mesh_1.GetGroups()
Sub_mesh_1 = Regular_1D.GetSubMesh()
Sub_mesh_2 = Regular_1D_1.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Edge_1_1, 'Edge_1')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Edge_4_1, 'Edge_4')
smesh.SetName(Edge_2_1, 'Edge_2')
smesh.SetName(Edge_3_1, 'Edge_3')

###
### HEXABLOCK component
###


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
