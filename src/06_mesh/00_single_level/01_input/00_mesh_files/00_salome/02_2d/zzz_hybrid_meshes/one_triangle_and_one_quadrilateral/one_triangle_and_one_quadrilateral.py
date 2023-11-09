#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/src/06_mesh/00_single_level/01_input/00_mesh_files/00_salome/02_2d')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("lx", 1)
notebook.set("ly", 1)
notebook.set("origin_x", 0)
notebook.set("origin_y", 0)
notebook.set("vertex_x", "origin_x+ lx")
notebook.set("vertex_y", "origin_y + ly")
notebook.set("square_below_delta_y", -2)
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
Vertex_1 = geompy.MakeVertex("origin_x", "origin_y", 0)
Vertex_2 = geompy.MakeVertex("vertex_x", "origin_y", 0)
Vertex_3 = geompy.MakeVertex("origin_x", "vertex_y", 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_3 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_1)
Face_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_3], 1)
[Edge_1,Edge_2,Edge_3] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)
[Edge_1, Edge_2, Edge_3] = geompy.GetExistingSubObjects(Face_1, False)
Vertex_4 = geompy.MakeVertex("origin_x", "square_below_delta_y", 0)
Vertex_5 = geompy.MakeVertex("vertex_x", "square_below_delta_y", 0)
Line_4 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_5)
Line_5 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_2)
Line_6 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_4)
Face_2 = geompy.MakeFaceWires([Line_1, Line_4, Line_5, Line_6], 1)
Compound_1 = geompy.MakeCompound([Face_1, Face_2])
[Face_3,Face_4] = geompy.ExtractShapes(Compound_1, geompy.ShapeType["FACE"], True)
[geomObj_1,geomObj_2,geomObj_3,geomObj_4,geomObj_5,geomObj_6,geomObj_7] = geompy.ExtractShapes(Compound_1, geompy.ShapeType["EDGE"], True)
[Edge_8,Edge_9,Edge_10] = geompy.ExtractShapes(Face_3, geompy.ShapeType["EDGE"], True)
[Edge_4,Edge_5,Edge_6,Edge_7] = geompy.ExtractShapes(Face_4, geompy.ShapeType["EDGE"], True)
Auto_group_for_Group_6_0 = geompy.CreateGroup(Compound_1, geompy.ShapeType["FACE"])
geompy.UnionList(Auto_group_for_Group_6_0, [Face_3, Face_4])
Auto_group_for_Group_1_0 = geompy.CreateGroup(Compound_1, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_1_0, [Edge_8, Edge_4])
Auto_group_for_Group_2 = geompy.CreateGroup(Compound_1, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_2, [Edge_10, Edge_5, Edge_7])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudyInFather( Face_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_1, Edge_3, 'Edge_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudy( Line_6, 'Line_6' )
geompy.addToStudy( Face_2, 'Face_2' )
geompy.addToStudy( Compound_1, 'Compound_1' )
geompy.addToStudyInFather( Compound_1, Face_3, 'Face_3' )
geompy.addToStudyInFather( Compound_1, Face_4, 'Face_4' )
geompy.addToStudyInFather( Face_4, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( Face_4, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_4, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_4, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Compound_1, Auto_group_for_Group_6_0, 'Auto_group_for_Group_6_0' )
geompy.addToStudyInFather( Face_3, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Face_3, Edge_9, 'Edge_9' )
geompy.addToStudyInFather( Face_3, Edge_10, 'Edge_10' )
geompy.addToStudyInFather( Compound_1, Auto_group_for_Group_1_0, 'Auto_group_for_Group_1_0' )
geompy.addToStudyInFather( Compound_1, Auto_group_for_Group_2, 'Auto_group_for_Group_2' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

aFilterManager = smesh.CreateFilterManager()
NETGEN_2D_Simple_Parameters_1 = smesh.CreateHypothesis('NETGEN_SimpleParameters_2D', 'NETGENEngine')
NETGEN_2D_Simple_Parameters_1.SetNumberOfSegments( 1 )
NETGEN_2D_Simple_Parameters_1.LengthFromEdges()
NETGEN_2D_Simple_Parameters_1.SetAllowQuadrangles( 0 )
NETGEN_1D_2D = smesh.CreateHypothesis('NETGEN_2D', 'NETGENEngine')
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,Face_1)
aCriteria.append(aCriterion)
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,'Face_1')
aCriteria.append(aCriterion)
aFilter0x4fd71d0 = aFilterManager.CreateFilter()
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,'Face_1')
aCriteria.append(aCriterion)
try:
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
Mesh_1 = smesh.Mesh(Compound_1,'Mesh_1')
NETGEN_2D_Simple_Parameters_2 = smesh.CreateHypothesis('NETGEN_SimpleParameters_2D', 'NETGENEngine')
NETGEN_2D_Simple_Parameters_2.SetNumberOfSegments( 1 )
NETGEN_2D_Simple_Parameters_2.LengthFromEdges()
NETGEN_2D_Simple_Parameters_2.SetAllowQuadrangles( 0 )
status = Mesh_1.AddHypothesis(NETGEN_1D_2D,Face_1)
status = Mesh_1.AddHypothesis(NETGEN_2D_Simple_Parameters_2,Face_1)
Regular_1D = Mesh_1.Segment(geom=Edge_4)
Number_of_Segments_1 = Regular_1D.NumberOfSegments(1)
Regular_1D_1 = Mesh_1.Segment(geom=Edge_6)
status = Mesh_1.AddHypothesis(Number_of_Segments_1,Edge_6)
Regular_1D_2 = Mesh_1.Segment(geom=Face_4)
Propagation_of_1D_Hyp = Regular_1D_2.Propagation()
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE,geom=Face_4)
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.EDGE,SMESH.FT_FreeBorders,SMESH.FT_Undefined,0,SMESH.FT_LogicalNOT)
aCriteria.append(aCriterion)
#smeshObj_1 = Mesh_1.GetMesh().UnionListOfGroups([ smeshObj_2, smeshObj_3 ], 'Group_6_0' ) ### smeshObj_2 has not been yet created
#Mesh_1.RemoveGroup( smeshObj_1 ) ### smeshObj_1 has not been yet created
Group_6_0 = Mesh_1.GroupOnGeom(Auto_group_for_Group_6_0,'Group_6_0',SMESH.FACE)
[ Group_6_0 ] = Mesh_1.GetGroups()
Group_1_0 = Mesh_1.GroupOnGeom(Auto_group_for_Group_1_0,'Group_1_0',SMESH.EDGE)
[ Group_6_0, Group_1_0 ] = Mesh_1.GetGroups()
Group_2_0 = Mesh_1.GroupOnGeom(Auto_group_for_Group_2,'Group_2',SMESH.EDGE)
[ Group_6_0, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Group_2_0.SetName( 'Group_2_0' )
[ Group_6_0, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
isDone = Mesh_1.Compute()
[ Group_6_0, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
coincident_nodes_on_part = Mesh_1.FindCoincidentNodesOnPart( [ Mesh_1 ], 1e-05, [], 0 )
Mesh_1.MergeNodes([[ 2, 5 ], [ 3, 4 ]], [], 0)
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.EDGE,SMESH.FT_FreeBorders,SMESH.FT_Undefined,0,SMESH.FT_LogicalNOT)
aCriteria.append(aCriterion)
isDone = Mesh_1.RemoveElements( [ 1, 5 ] )
[ Group_6_0, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0)
[ Group_6_0, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertFromQuadratic()
[ Group_6_0, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
[ Group_6_0, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0)
[ Group_6_0, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED( r'/home/gbornia/software/femus/src/06_mesh/00_single_level/01_input/00_mesh_files/00_salome/02_2d/one_tri6_one_quad8_groups_b_v.med', 0, 41, 1, Mesh_1, 0, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
[ Group_6_0, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED( r'/home/gbornia/software/femus/src/06_mesh/00_single_level/01_input/00_mesh_files/00_salome/02_2d/one_tri7_one_quad9_groups_b_v.med', 0, 41, 1, Mesh_1, 0, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
Sub_mesh_1 = Mesh_1.GetSubMesh( Face_3, 'Sub-mesh_1' )
Sub_mesh_2 = Regular_1D.GetSubMesh()
Sub_mesh_3 = Regular_1D_1.GetSubMesh()
Sub_mesh_4 = Regular_1D_2.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(NETGEN_1D_2D, 'NETGEN 1D-2D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(NETGEN_2D_Simple_Parameters_1, 'NETGEN 2D Simple Parameters_1')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(NETGEN_2D_Simple_Parameters_2, 'NETGEN 2D Simple Parameters_2')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Group_6_0, 'Group_6_0')
smesh.SetName(Group_2_0, 'Group_2_0')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
