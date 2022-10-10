#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.8.0 with dump python functionality
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
notebook.set("geom_l_x", 1)
notebook.set("geom_l_y", 1)
notebook.set("geom_x_b", 0)
notebook.set("geom_x_e", "geom_x_b + geom_l_x")
notebook.set("geom_y_b", 0)
notebook.set("geom_y_e", "geom_y_b + geom_l_y")
notebook.set("geom_delta_x_from_x_b", 0.75)
notebook.set("geom_delta_y_from_y_b", 0.8)
notebook.set("geom_indent_l_x", "geom_l_x - geom_delta_x_from_x_b")
notebook.set("geom_indent_l_y", "geom_l_y - geom_delta_y_from_y_b")
notebook.set("mesh_indent_l_x", 1)
notebook.set("mesh_indent_l_y", 1)
notebook.set("mesh_delta_x_from_x_b", 1)
notebook.set("mesh_delta_y_from_y_b", 1)
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
Vertex_1 = geompy.MakeVertex("geom_x_b", "geom_y_b", 0)
Vertex_2 = geompy.MakeVertex("geom_x_b", "geom_y_e", 0)
Vertex_3 = geompy.MakeVertex("geom_x_e", "geom_y_e", 0)
Vertex_4 = geompy.MakeVertex("geom_x_e", "geom_y_b", 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_3 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_4)
Line_4 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_1)
Face_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_3, Line_4], 1)
Vertex_5 = geompy.MakeVertex("geom_delta_x_from_x_b", "geom_delta_y_from_y_b", 0)
Vertex_6 = geompy.MakeVertex("geom_delta_x_from_x_b", "geom_y_e", 0)
Vertex_7 = geompy.MakeVertex("geom_x_e", "geom_delta_y_from_y_b", 0)
Vertex_8 = geompy.MakeVertex("geom_x_e", "geom_y_e", 0)
Edge_1 = geompy.MakeEdge(Vertex_5, Vertex_7)
Edge_2 = geompy.MakeEdge(Vertex_7, Vertex_8)
Edge_3 = geompy.MakeEdge(Vertex_8, Vertex_6)
Edge_4 = geompy.MakeEdge(Vertex_6, Vertex_5)
Face_2 = geompy.MakeFaceWires([Edge_1, Edge_2, Edge_3, Edge_4], 1)
geomObj_1 = geompy.MakeCutList(Face_1, [Face_2], True)
[geomObj_2,geomObj_3,geomObj_4,geomObj_5,geomObj_6,geomObj_7] = geompy.ExtractShapes(geomObj_1, geompy.ShapeType["EDGE"], True)
[geomObj_2, geomObj_3, geomObj_4, geomObj_5, geomObj_6, geomObj_7] = geompy.GetExistingSubObjects(geomObj_1, False)
Cut_1 = geompy.MakeCutList(Face_1, [Face_2], True)
Vertex_9 = geompy.MakeVertex("geom_delta_x_from_x_b", "geom_y_b", 0)
Vertex_10 = geompy.MakeVertex("geom_x_b", "geom_delta_y_from_y_b", 0)
Line_5 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_9)
Line_6 = geompy.MakeLineTwoPnt(Vertex_10, Vertex_5)
Partition_1 = geompy.MakePartition([Cut_1], [Line_5, Line_6], [], [], geompy.ShapeType["FACE"], 0, [], 0)
[Edge_5,Edge_6,Edge_7,Edge_8,Edge_9,Edge_10,Edge_11,Edge_12,Edge_13,Edge_14] = geompy.ExtractShapes(Partition_1, geompy.ShapeType["EDGE"], True)
Auto_group_for_Group_1_0 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_1_0, [Edge_5, Edge_6])
Auto_group_for_Group_3_0 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_3_0, [Edge_7, Edge_12])
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
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Vertex_7, 'Vertex_7' )
geompy.addToStudy( Vertex_8, 'Vertex_8' )
geompy.addToStudy( Edge_1, 'Edge_1' )
geompy.addToStudy( Edge_2, 'Edge_2' )
geompy.addToStudy( Edge_3, 'Edge_3' )
geompy.addToStudy( Edge_4, 'Edge_4' )
geompy.addToStudy( Face_2, 'Face_2' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudy( Vertex_9, 'Vertex_9' )
geompy.addToStudy( Vertex_10, 'Vertex_10' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudy( Line_6, 'Line_6' )
geompy.addToStudy( Partition_1, 'Partition_1' )
geompy.addToStudyInFather( Partition_1, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Partition_1, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Partition_1, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Partition_1, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Partition_1, Edge_9, 'Edge_9' )
geompy.addToStudyInFather( Partition_1, Edge_10, 'Edge_10' )
geompy.addToStudyInFather( Partition_1, Edge_11, 'Edge_11' )
geompy.addToStudyInFather( Partition_1, Edge_12, 'Edge_12' )
geompy.addToStudyInFather( Partition_1, Edge_13, 'Edge_13' )
geompy.addToStudyInFather( Partition_1, Edge_14, 'Edge_14' )
geompy.addToStudyInFather( Partition_1, Auto_group_for_Group_1_0, 'Auto_group_for_Group_1_0' )
geompy.addToStudyInFather( Partition_1, Auto_group_for_Group_3_0, 'Auto_group_for_Group_3_0' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Number_of_Segments_1 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_1.SetNumberOfSegments( "mesh_indent_l_y" )
Number_of_Segments_1.SetReversedEdges( [] )
Number_of_Segments_2 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_2.SetNumberOfSegments( "mesh_indent_l_x" )
Number_of_Segments_2.SetReversedEdges( [] )
Number_of_Segments_3 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_3.SetNumberOfSegments( "mesh_delta_y_from_y_b" )
Number_of_Segments_3.SetReversedEdges( [] )
Number_of_Segments_4 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_4.SetNumberOfSegments( "mesh_delta_x_from_x_b" )
Number_of_Segments_4.SetReversedEdges( [] )
Number_of_Segments_5 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_5.SetNumberOfSegments( 15 )
Number_of_Segments_5.SetReversedEdges( [] )
Mesh_1 = smesh.Mesh(Partition_1)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D = Mesh_1.Segment(geom=Edge_5)
Number_of_Segments_6 = Regular_1D.NumberOfSegments("mesh_delta_y_from_y_b",None,[])
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Number_of_Segments_7 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_7.SetNumberOfSegments( "mesh_indent_l_x" )
Number_of_Segments_8 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_8.SetNumberOfSegments( "mesh_delta_x_from_x_b" )
Regular_1D_1 = Mesh_1.Segment(geom=Edge_12)
Number_of_Segments_9 = Regular_1D_1.NumberOfSegments("mesh_indent_l_x",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_12)
Propagation_of_1D_Hyp_1 = smesh.CreateHypothesis('Propagation')
Regular_1D_2 = Mesh_1.Segment(geom=Edge_6)
Number_of_Segments_10 = Regular_1D_2.NumberOfSegments("mesh_indent_l_y",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_6)
Regular_1D_3 = Mesh_1.Segment(geom=Edge_7)
Number_of_Segments_11 = Regular_1D_3.NumberOfSegments("mesh_delta_x_from_x_b",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_7)
Regular_1D_4 = Mesh_1.Segment()
isDone = Mesh_1.Compute()
Group_1_0 = Mesh_1.GroupOnGeom(Auto_group_for_Group_1_0,'Group_1_0',SMESH.EDGE)
[ Group_1_0 ] = Mesh_1.GetGroups()
Group_2_0 = Mesh_1.GroupOnGeom(Edge_14,'Group_2_0',SMESH.EDGE)
[ Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Group_3_0 = Mesh_1.GroupOnGeom(Auto_group_for_Group_3_0,'Group_3_0',SMESH.EDGE)
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
Group_4_0 = Mesh_1.GroupOnGeom(Edge_9,'Group_4_0',SMESH.EDGE)
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0 ] = Mesh_1.GetGroups()
Group_5_0 = Mesh_1.GroupOnGeom(Edge_11,'Group_5_0',SMESH.EDGE)
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0, Group_5_0 ] = Mesh_1.GetGroups()
Group_6_0 = Mesh_1.GroupOnGeom(Edge_13,'Group_6_0',SMESH.EDGE)
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0, Group_5_0, Group_6_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0, Group_5_0, Group_6_0 ] = Mesh_1.GetGroups()
Sub_mesh_1 = Regular_1D.GetSubMesh()
Sub_mesh_2 = Regular_1D_1.GetSubMesh()
Sub_mesh_3 = Regular_1D_2.GetSubMesh()
Sub_mesh_4 = Regular_1D_3.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Number_of_Segments_10, 'Number of Segments_10')
smesh.SetName(Number_of_Segments_11, 'Number of Segments_11')
smesh.SetName(Number_of_Segments_9, 'Number of Segments_9')
smesh.SetName(Propagation_of_1D_Hyp_1, 'Propagation of 1D Hyp. on Opposite Edges_2')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Group_5_0, 'Group_5_0')
smesh.SetName(Group_6_0, 'Group_6_0')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_6, 'Number of Segments_6')
smesh.SetName(Number_of_Segments_5, 'Number of Segments_5')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Number_of_Segments_8, 'Number of Segments_8')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Number_of_Segments_7, 'Number of Segments_7')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
