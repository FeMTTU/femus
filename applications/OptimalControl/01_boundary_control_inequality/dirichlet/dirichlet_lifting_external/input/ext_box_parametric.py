#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.2.1 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/applications/OptimalControl/boundary_control_inequality/dirichlet/dirichlet_lifting_external/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("l_x", 1)
notebook.set("l_y", 1)
notebook.set("n_x", 5)
notebook.set("n_y", 4)
notebook.set("x_b", 0)
notebook.set("x_e", "x_b + l_x")
notebook.set("y_b", 0)
notebook.set("y_e", "y_b + l_y")
notebook.set("z_b", 0)
notebook.set("z_e", 0)
notebook.set("ext_l_x", 0.25)
notebook.set("ext_l_y", 0.5)
notebook.set("ext_x_b", "x_e")
notebook.set("ext_x_e", "ext_x_b + ext_l_x")
notebook.set("ext_y_b", 0.25)
notebook.set("ext_y_e", "ext_y_b + ext_l_y")
notebook.set("ext_n_x", 4)
notebook.set("ext_n_y", 2)
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
aFilterManager = smesh.CreateFilterManager()
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
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Group_4_0 = Mesh_1.GroupOnGeom(Edge_1,'Edge_1',SMESH.EDGE)
Group_4_0.SetName( 'Group_4_0' )
Group_1_0 = Mesh_1.GroupOnGeom(Edge_2,'Group_1',SMESH.EDGE)
Group_1_0.SetName( 'Group_1_0' )
Group_2_0 = Mesh_1.GroupOnGeom(Edge_4,'Group_2_0',SMESH.EDGE)
Group_3_0 = Mesh_1.GroupOnGeom(Edge_3,'Group_3_0',SMESH.EDGE)
Group_12_0 = Mesh_1.GroupOnGeom(Face_1,'Group_12_0',SMESH.FACE)
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/gbornia/software/femus/applications/OptimalControl/boundary_control_inequality/dirichlet/dirichlet_boundary/input/square_parametric.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Mesh_2 = smesh.Mesh(Face_2)
Regular_1D_3 = Mesh_2.Segment()
Quadrangle_2D_1 = Mesh_2.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_4 = Mesh_2.Segment(geom=Edge_5)
Number_of_Segments_3 = Regular_1D_4.NumberOfSegments("ext_n_y",None,[])
status = Mesh_2.AddHypothesis(Propagation_of_1D_Hyp,Edge_5)
Regular_1D_5 = Mesh_2.Segment(geom=Edge_6)
Number_of_Segments_4 = Regular_1D_5.NumberOfSegments("ext_n_x",None,[])
status = Mesh_2.AddHypothesis(Propagation_of_1D_Hyp,Edge_6)
isDone = Mesh_2.Compute()
Mesh_2.ConvertToQuadratic(0, Mesh_2,True)
[ Group_4_0, Group_1_0, Group_2_0, Group_3_0, Group_12_0 ] = Mesh_1.GetGroups()
Mesh_3 = smesh.Concatenate([ Mesh_1.GetMesh(), Mesh_2.GetMesh() ], 1, 1, 1e-05)
aFilterLibrary0x9a6fa80 = aFilterManager.LoadLibrary('/home/gbornia/FilterLib.xml')
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.EDGE,SMESH.FT_FreeBorders,SMESH.FT_Undefined,0,SMESH.FT_LogicalNOT,SMESH.FT_Undefined,0); aCriterion.Precision = 0
aCriteria.append(aCriterion)
aFilterLibrary0x9a6fa80.Add('EdgeFilter_1',aFilter0xd437850)
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.EDGE,SMESH.FT_FreeBorders,SMESH.FT_Undefined,0,SMESH.FT_LogicalNOT)
aCriteria.append(aCriterion)
aFilterLibrary0x9a6fa80.Copy('EdgeFilter_1')
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.EDGE,SMESH.FT_FreeBorders,SMESH.FT_Undefined,0,SMESH.FT_LogicalNOT)
aCriteria.append(aCriterion)
isDone = Mesh_3.RemoveElements( [ 7, 8 ] )
#Mesh_3.RemoveGroup( smeshObj_1 ) ### smeshObj_1 has not been yet created
#Mesh_3.RemoveGroup( smeshObj_2 ) ### smeshObj_2 has not been yet created
#Mesh_3.RemoveGroup( smeshObj_3 ) ### smeshObj_3 has not been yet created
#Mesh_3.RemoveGroup( smeshObj_4 ) ### smeshObj_4 has not been yet created
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,'Face_1')
aCriteria.append(aCriterion)
aFilter_4 = smesh.GetFilterFromCriteria(aCriteria)
aFilter_4.SetMesh(Mesh_3.GetMesh())
Group_12_0_1 = Mesh_3.GroupOnFilter( SMESH.FACE, 'Group_12_0', aFilter_4 )
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,'Face_2')
aCriteria.append(aCriterion)
aFilter_5 = smesh.GetFilterFromCriteria(aCriteria)
aFilter_5.SetMesh(Mesh_3.GetMesh())
Group_13_0 = Mesh_3.GroupOnFilter( SMESH.FACE, 'Group_13_0', aFilter_5 )
Group_1_0_1 = Mesh_3.CreateEmptyGroup( SMESH.EDGE, 'Group_1_0' )
nbAdd = Group_1_0_1.AddFrom( Mesh_3.GetMesh() )
Sub_mesh_1 = Regular_1D.GetSubMesh()
Sub_mesh_2 = Regular_1D_1.GetSubMesh()
Sub_mesh_3 = Regular_1D_4.GetSubMesh()
Sub_mesh_4 = Regular_1D_5.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Group_13_0, 'Group_13_0')
smesh.SetName(Group_12_0_1, 'Group_12_0')
smesh.SetName(Group_1_0_1, 'Group_1_0')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Group_12_0, 'Group_12_0')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')
smesh.SetName(Mesh_3.GetMesh(), 'Mesh_3')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Propagation_of_1D_Hyp_1, 'Propagation of 1D Hyp. on Opposite Edges_2')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_3')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
