#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.3.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("l_x", 1)
notebook.set("l_y", 1)
notebook.set("leg_distance_from_x_axis", 0.5)
notebook.set("leg_distance_from_y_axis", 0.5)
notebook.set("mesh_n_x_core", 1)
notebook.set("mesh_n_y_core", 1)
notebook.set("mesh_n_x_leg", 1)
notebook.set("mesh_n_y_leg", 1)
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
Vertex_1 = geompy.MakeVertex(0, 0, 0)
Vertex_3 = geompy.MakeVertex("l_x", 0, 0)
Vertex_2 = geompy.MakeVertex(0, "l_y", 0)
Vertex_4 = geompy.MakeVertex("leg_distance_from_x_axis", "leg_distance_from_y_axis", 0)
Vertex_5 = geompy.MakeVertex("leg_distance_from_x_axis", "l_y", 0)
Vertex_6 = geompy.MakeVertex("l_x", "leg_distance_from_y_axis", 0)
Vertex_7 = geompy.MakeVertex("leg_distance_from_x_axis", 0, 0)
Vertex_8 = geompy.MakeVertex(0, "leg_distance_from_y_axis", 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_7)
Line_1_vertex_3 = geompy.GetSubShape(Line_1, [3])
Line_2 = geompy.MakeLineTwoPnt(Vertex_4, Line_1_vertex_3)
Line_3 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_8)
Line_1_vertex_2 = geompy.GetSubShape(Line_1, [2])
Line_4 = geompy.MakeLineTwoPnt(Vertex_8, Line_1_vertex_2)
Line_3_vertex_2 = geompy.GetSubShape(Line_3, [2])
Line_5 = geompy.MakeLineTwoPnt(Line_3_vertex_2, Vertex_5)
Line_6 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_2)
Line_3_vertex_3 = geompy.GetSubShape(Line_3, [3])
Line_7 = geompy.MakeLineTwoPnt(Vertex_2, Line_3_vertex_3)
Line_8 = geompy.MakeLineTwoPnt(Line_1_vertex_3, Vertex_3)
Line_9 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_6)
Line_10 = geompy.MakeLineTwoPnt(Vertex_6, Line_3_vertex_2)
Face_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_3, Line_4], 1)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)
Face_2 = geompy.MakeFaceWires([Line_3, Line_5, Line_6, Line_7], 1)
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Face_2, geompy.ShapeType["EDGE"], True)
Face_3 = geompy.MakeFaceWires([Line_2, Line_8, Line_9, Line_10], 1)
[Edge_9,Edge_10,Edge_11,Edge_12] = geompy.ExtractShapes(Face_3, geompy.ShapeType["EDGE"], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Vertex_7, 'Vertex_7' )
geompy.addToStudy( Vertex_8, 'Vertex_8' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudyInFather( Line_1, Line_1_vertex_3, 'Line_1:vertex_3' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudyInFather( Line_1, Line_1_vertex_2, 'Line_1:vertex_2' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudyInFather( Line_3, Line_3_vertex_2, 'Line_3:vertex_2' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudy( Line_6, 'Line_6' )
geompy.addToStudyInFather( Line_3, Line_3_vertex_3, 'Line_3:vertex_3' )
geompy.addToStudy( Line_7, 'Line_7' )
geompy.addToStudy( Line_8, 'Line_8' )
geompy.addToStudy( Line_9, 'Line_9' )
geompy.addToStudy( Line_10, 'Line_10' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudyInFather( Face_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_1, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Face_1, Edge_4, 'Edge_4' )
geompy.addToStudy( Face_2, 'Face_2' )
geompy.addToStudyInFather( Face_2, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_2, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_2, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Face_2, Edge_8, 'Edge_8' )
geompy.addToStudy( Face_3, 'Face_3' )
geompy.addToStudyInFather( Face_3, Edge_9, 'Edge_9' )
geompy.addToStudyInFather( Face_3, Edge_10, 'Edge_10' )
geompy.addToStudyInFather( Face_3, Edge_11, 'Edge_11' )
geompy.addToStudyInFather( Face_3, Edge_12, 'Edge_12' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Face_1)
Mesh_2 = smesh.Mesh(Face_2)
Mesh_3 = smesh.Mesh(Face_3)
Regular_1D = Mesh_1.Segment(geom=Edge_2)
Number_of_Segments_1 = Regular_1D.NumberOfSegments("mesh_n_x_core",None,[])
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Regular_1D_1 = Mesh_1.Segment(geom=Edge_1)
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_1)
Number_of_Segments_2 = Regular_1D_1.NumberOfSegments("mesh_n_y_core",None,[])
Regular_1D_3 = Mesh_2.Segment(geom=Edge_5)
Number_of_Segments_3 = Regular_1D_3.NumberOfSegments("mesh_n_y_leg",None,[])
status = Mesh_2.AddHypothesis(Propagation_of_1D_Hyp,Edge_5)
Regular_1D_4 = Mesh_2.Segment(geom=Edge_6)
Number_of_Segments_4 = Regular_1D_4.NumberOfSegments("mesh_n_x_core",None,[])
status = Mesh_2.AddHypothesis(Propagation_of_1D_Hyp,Edge_6)
Regular_1D_6 = Mesh_3.Segment(geom=Edge_9)
Number_of_Segments_5 = Regular_1D_6.NumberOfSegments("mesh_n_y_core",None,[])
status = Mesh_3.AddHypothesis(Propagation_of_1D_Hyp,Edge_9)
Regular_1D_7 = Mesh_3.Segment(geom=Edge_10)
Number_of_Segments_6 = Regular_1D_7.NumberOfSegments("mesh_n_x_leg",None,[])
status = Mesh_3.AddHypothesis(Propagation_of_1D_Hyp,Edge_10)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_2 = Mesh_1.Segment()
isDone = Mesh_1.Compute()
Regular_1D_5 = Mesh_2.Segment()
Quadrangle_2D_1 = Mesh_2.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_2.Compute()
Regular_1D_8 = Mesh_3.Segment()
Quadrangle_2D_2 = Mesh_3.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_3.Compute()
Mesh_3.ConvertToQuadratic(0, Mesh_3,True)
Mesh_2.ConvertToQuadratic(0, Mesh_2,True)
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Mesh_4 = smesh.Concatenate( [ Mesh_1.GetMesh(), Mesh_2.GetMesh(), Mesh_3.GetMesh() ], 1, 1, 1e-05, False )
isDone = Mesh_4.RemoveElements( [ 2, 3 ] )
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.EDGE,SMESH.FT_FreeBorders,SMESH.FT_Undefined,0,SMESH.FT_LogicalNOT)
aCriteria.append(aCriterion)
Group_3_0 = Mesh_4.CreateEmptyGroup( SMESH.EDGE, 'Group_1_0' )
nbAdd = Group_3_0.Add( [ 1, 8 ] )
Group_4_0 = Mesh_4.CreateEmptyGroup( SMESH.EDGE, 'Group_2_0' )
nbAdd = Group_4_0.Add( [ 5 ] )
Group_1_0 = Mesh_4.CreateEmptyGroup( SMESH.EDGE, 'Group_3_0' )
nbAdd = Group_1_0.Add( [ 2, 4 ] )
Group_2_0 = Mesh_4.CreateEmptyGroup( SMESH.EDGE, 'Group_4_0' )
nbAdd = Group_2_0.Add( [ 9 ] )
Group_3_0.SetName( 'Group_13_0' )
Group_4_0.SetName( 'Group_24_0' )
Group_1_0.SetName( 'Group_1_0' )
Group_2_0.SetName( 'Group_2_0' )
Group_3_0.SetName( 'Group_3_0' )
Group_4_0.SetName( 'Group_4_0' )
Group_5_0 = Mesh_4.CreateEmptyGroup( SMESH.EDGE, 'Group_5_0' )
nbAdd = Group_5_0.Add( [ 6 ] )
Group_6_0 = Mesh_4.CreateEmptyGroup( SMESH.EDGE, 'Group_6_0' )
nbAdd = Group_6_0.Add( [ 10 ] )
Group_7_0 = Mesh_4.CreateEmptyGroup( SMESH.FACE, 'Group_7_0' )
nbAdd = Group_7_0.Add( [ 3, 7, 11 ] )
smesh.SetName(Mesh_4, 'Mesh_4')
try:
  Mesh_4.ExportMED(r'/home/gbornia/Lshape.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Sub_mesh_1 = Regular_1D.GetSubMesh()
Sub_mesh_2 = Regular_1D_1.GetSubMesh()
Sub_mesh_3 = Regular_1D_3.GetSubMesh()
Sub_mesh_4 = Regular_1D_4.GetSubMesh()
Sub_mesh_5 = Regular_1D_6.GetSubMesh()
Sub_mesh_6 = Regular_1D_7.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Group_7_0, 'Group_7_0')
smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
smesh.SetName(Sub_mesh_6, 'Sub-mesh_6')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Group_5_0, 'Group_5_0')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')
smesh.SetName(Mesh_3.GetMesh(), 'Mesh_3')
smesh.SetName(Group_6_0, 'Group_6_0')
smesh.SetName(Mesh_4.GetMesh(), 'Mesh_4')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_6, 'Number of Segments_6')
smesh.SetName(Number_of_Segments_5, 'Number of Segments_5')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
