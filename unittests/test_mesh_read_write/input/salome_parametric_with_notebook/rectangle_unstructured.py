#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.8.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/unittests/test_mesh_read_write/input/salome_parametric_with_notebook')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("l_x", 1)
notebook.set("l_y", 1)
notebook.set("n_x", 2)
notebook.set("n_y", 2)
notebook.set("x_b", 0)
notebook.set("x_e", "x_b + l_x")
notebook.set("y_b", 0)
notebook.set("y_e", "y_b + l_y")
notebook.set("z_b", 0)
notebook.set("z_e", 0)
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
O_1 = geompy.MakeVertex(0, 0, 0)
OX_1 = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY_1 = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ_1 = geompy.MakeVectorDXDYDZ(0, 0, 1)
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
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( O_1, 'O' )
geompy.addToStudy( OX_1, 'OX' )
geompy.addToStudy( OY_1, 'OY' )
geompy.addToStudy( OZ_1, 'OZ' )
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

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Face_1)
Regular_1D = Mesh_1.Segment(geom=Edge_2)
Number_of_Segments_1 = Regular_1D.NumberOfSegments("n_x",None,[])
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Propagation_of_1D_Hyp_1 = smesh.CreateHypothesis('Propagation')
Regular_1D_1 = Mesh_1.Segment(geom=Edge_1)
Number_of_Segments_2 = Regular_1D_1.NumberOfSegments("n_y",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_1)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_2 = Mesh_1.Segment()
Group_1_0 = Mesh_1.GroupOnGeom(Edge_1,'Edge_1',SMESH.EDGE)
Group_1_0.SetName( 'Group_4_0' )
Group_3_0 = Mesh_1.GroupOnGeom(Edge_2,'Group_1',SMESH.EDGE)
Group_3_0.SetName( 'Group_1_0' )
Group_2_0 = Mesh_1.GroupOnGeom(Edge_4,'Group_2_0',SMESH.EDGE)
Group_4_0 = Mesh_1.GroupOnGeom(Edge_3,'Group_3_0',SMESH.EDGE)
Group_12_0 = Mesh_1.GroupOnGeom(Face_1,'Group_12_0',SMESH.FACE)
[ Group_1_0, Group_3_0, Group_2_0, Group_4_0, Group_12_0 ] = Mesh_1.GetGroups()
Group_1_0.SetName( 'Group_1_0_n' )
Group_2_0.SetName( 'Group_3_0_n' )
Group_4_0.SetName( 'Group_4_0' )
Group_2_0.SetName( 'Group_2_0_n' )
Group_3_0.SetName( 'Group_3_0' )
Group_2_0.SetName( 'Group_2_0' )
Group_1_0.SetName( 'Group_1_0' )
status = Mesh_1.RemoveHypothesis(Number_of_Segments_1,Edge_2)
Number_of_Segments_3 = Regular_1D.NumberOfSegments("n_x",1.3,[])
[ Group_1_0, Group_3_0, Group_2_0, Group_4_0, Group_12_0 ] = Mesh_1.GetGroups()
status = Mesh_1.RemoveHypothesis(Number_of_Segments_3,Edge_2)
Number_of_Segments_4 = Regular_1D.NumberOfSegments("n_x",1.8,[])
[ Group_1_0, Group_3_0, Group_2_0, Group_4_0, Group_12_0 ] = Mesh_1.GetGroups()
Propagation_of_Node = smesh.CreateHypothesis('PropagOfDistribution')
Propagation_of_1D_Hyp_2 = smesh.CreateHypothesis('Propagation')
status = Mesh_1.RemoveHypothesis(Propagation_of_1D_Hyp,Edge_2)
Propagation_of_Node_1 = Regular_1D.PropagationOfDistribution()
[ Group_1_0, Group_3_0, Group_2_0, Group_4_0, Group_12_0 ] = Mesh_1.GetGroups()
Regular_1D_3 = Mesh_1.Segment(geom=Edge_3)
Number_of_Segments_5 = Regular_1D_3.NumberOfSegments(2)
[ Group_1_0, Group_3_0, Group_2_0, Group_4_0, Group_12_0 ] = Mesh_1.GetGroups()
isDone = Mesh_1.Compute()
[ Group_1_0, Group_3_0, Group_2_0, Group_4_0, Group_12_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
[ Group_1_0, Group_3_0, Group_2_0, Group_4_0, Group_12_0 ] = Mesh_1.GetGroups()
Sub_mesh_1 = Regular_1D.GetSubMesh()
Sub_mesh_2 = Regular_1D_1.GetSubMesh()
Sub_mesh_3 = Regular_1D_3.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Propagation_of_1D_Hyp_2, 'Propagation of 1D Hyp. on Opposite Edges_3')
smesh.SetName(Number_of_Segments_5, 'Number of Segments_5')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Group_12_0, 'Group_12_0')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Propagation_of_1D_Hyp_1, 'Propagation of 1D Hyp. on Opposite Edges_2')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Propagation_of_Node, 'Propagation of Node Distribution on Opposite Edges_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Propagation_of_Node_1, 'Propagation of Node Distribution on Opposite Edges_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
