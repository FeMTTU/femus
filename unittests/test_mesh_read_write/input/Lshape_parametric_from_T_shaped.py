#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.6.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/unittests/test_mesh_read_write/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("l_x", 0.5)
notebook.set("l_y", 1)
notebook.set("n_x", 1)
notebook.set("n_y", 2)
notebook.set("x_b", 0)
notebook.set("x_e", "x_b + l_x")
notebook.set("y_b", 0)
notebook.set("y_e", "y_b + l_y")
notebook.set("z_b", 0)
notebook.set("z_e", 0)
notebook.set("ext_l_x", 0.5)
notebook.set("ext_l_y", 0.5)
notebook.set("ext_x_b", "x_e")
notebook.set("ext_x_e", "ext_x_b + ext_l_x")
notebook.set("ext_y_b", 0)
notebook.set("ext_y_e", "ext_y_b + ext_l_y")
notebook.set("ext_n_x", 1)
notebook.set("ext_n_y", 1)
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
Vertex_1 = geompy.MakeVertex("x_b", "y_b", "z_b")
Vertex_2 = geompy.MakeVertex("x_b", "y_e", 0)
Vertex_3 = geompy.MakeVertex("x_e", "y_b", 0)
Vertex_4 = geompy.MakeVertex("x_e", "y_e", 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_3)
Line_2 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_4)
Line_3 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_2)
Line_4 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_1)
Face_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_3, Line_4], 1)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)
Vertex_5 = geompy.MakeVertex("ext_x_b", "ext_y_b", 0)
Vertex_6 = geompy.MakeVertex("ext_x_e", "ext_y_b", 0)
Vertex_7 = geompy.MakeVertex("ext_x_e", "ext_y_e", 0)
Vertex_8 = geompy.MakeVertex("ext_x_b", "ext_y_e", 0)
Line_5 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_6)
Line_6 = geompy.MakeLineTwoPnt(Vertex_6, Vertex_7)
Line_7 = geompy.MakeLineTwoPnt(Vertex_7, Vertex_8)
Line_5_vertex_2 = geompy.GetSubShape(Line_5, [2])
Line_8 = geompy.MakeLineTwoPnt(Vertex_8, Line_5_vertex_2)
Face_2 = geompy.MakeFaceWires([Line_5, Line_6, Line_7, Line_8], 1)
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Face_2, geompy.ShapeType["EDGE"], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudyInFather( Face_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_1, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Face_1, Edge_4, 'Edge_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Vertex_7, 'Vertex_7' )
geompy.addToStudy( Vertex_8, 'Vertex_8' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudy( Line_6, 'Line_6' )
geompy.addToStudy( Line_7, 'Line_7' )
geompy.addToStudyInFather( Line_5, Line_5_vertex_2, 'Line_5:vertex_2' )
geompy.addToStudy( Line_8, 'Line_8' )
geompy.addToStudy( Face_2, 'Face_2' )
geompy.addToStudyInFather( Face_2, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_2, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_2, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Face_2, Edge_8, 'Edge_8' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

aFilterManager = smesh.CreateFilterManager()
Mesh_1 = smesh.Mesh(Face_1)
Regular_1D_1 = Mesh_1.Segment(geom=Edge_2)
NumberOfSegments_n_x = Regular_1D_1.NumberOfSegments("n_x",None,[])
Propagation = Regular_1D_1.Propagation()
Propagation_3 = smesh.CreateHypothesis('Propagation')
Regular_1D_1_1 = Mesh_1.Segment(geom=Edge_1)
NumberOfSegments_n_y = Regular_1D_1_1.NumberOfSegments("n_y",None,[])
status = Mesh_1.AddHypothesis(Propagation,Edge_1)
Quadrangle_2D_2 = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_1_2 = Mesh_1.Segment()
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Group_4_0 = Mesh_1.GroupOnGeom(Edge_1,'Edge_1',SMESH.EDGE)
Group_4_0.SetName( 'Group_4_0' )
Group_1_0 = Mesh_1.GroupOnGeom(Edge_2,'Group_1',SMESH.EDGE)
Group_1_0.SetName( 'Group_1_0' )
Group_2_0 = Mesh_1.GroupOnGeom(Edge_4,'Group_2_0',SMESH.EDGE)
Group_3_0 = Mesh_1.GroupOnGeom(Edge_3,'Group_3_0',SMESH.EDGE)
Group_12_0 = Mesh_1.GroupOnGeom(Face_1,'Group_12_0',SMESH.FACE)
Face_2_1 = smesh.Mesh(Face_2)
Regular_1D_1_3 = Face_2_1.Segment()
Quadrangle_2D_2_1 = Face_2_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_1_4 = Face_2_1.Segment(geom=Edge_5)
NumberOfSegments_ext = Regular_1D_1_4.NumberOfSegments("ext_n_y",None,[])
status = Face_2_1.AddHypothesis(Propagation,Edge_5)
Regular_1D_1_5 = Face_2_1.Segment(geom=Edge_6)
NumberOfSegments_ext_1 = Regular_1D_1_5.NumberOfSegments("ext_n_x",None,[])
status = Face_2_1.AddHypothesis(Propagation,Edge_6)
isDone = Face_2_1.Compute()
Face_2_1.ConvertToQuadratic(0, Face_2_1,True)
[ Group_4_0, Group_1_0, Group_2_0, Group_3_0, Group_12_0 ] = Mesh_1.GetGroups()
Mesh_5 = smesh.Concatenate( [ Mesh_1.GetMesh(), Face_2_1.GetMesh() ], 1, 1, 1e-05, False )
[ Group_4_0_1, Group_1_0_1, Group_2_0_1, Group_3_0_1, Group_12_0_1 ] = Mesh_5.GetGroups()
aFilterLibrary0x2971740 = aFilterManager.LoadLibrary('/home/gbornia/FilterLib.xml')
Regular_1D = Regular_1D_1.GetSubMesh()
Regular_1D_2 = Regular_1D_1_1.GetSubMesh()
Regular_1D_3 = Regular_1D_1_4.GetSubMesh()
Regular_1D_4 = Regular_1D_1_5.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Group_12_0_1, 'Group_12_0')
smesh.SetName(Group_3_0_1, 'Group_3_0')
smesh.SetName(Group_1_0_1, 'Group_1_0')
smesh.SetName(Group_2_0_1, 'Group_2_0')
smesh.SetName(Group_4_0_1, 'Group_4_0')
smesh.SetName(Regular_1D_1.GetAlgorithm(), 'Regular_1D_1')
smesh.SetName(Quadrangle_2D_2.GetAlgorithm(), 'Quadrangle_2D_2')
smesh.SetName(Group_12_0, 'Group_12_0')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Face_2_1.GetMesh(), 'Face_2')
smesh.SetName(Mesh_5.GetMesh(), 'Mesh_5')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Propagation_3, 'Propagation_3')
smesh.SetName(Propagation, 'Propagation')
smesh.SetName(NumberOfSegments_n_x, 'NumberOfSegments=n_x,[],0:1:1:13')
smesh.SetName(NumberOfSegments_ext_1, 'NumberOfSegments=ext_n_x,[],0:1:1:22')
smesh.SetName(NumberOfSegments_ext, 'NumberOfSegments=ext_n_y,[],0:1:1:22')
smesh.SetName(NumberOfSegments_n_y, 'NumberOfSegments=n_y,[],0:1:1:13')
smesh.SetName(Regular_1D_3, 'Regular_1D')
smesh.SetName(Regular_1D_2, 'Regular_1D')
smesh.SetName(Regular_1D_4, 'Regular_1D')
smesh.SetName(Regular_1D, 'Regular_1D')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
