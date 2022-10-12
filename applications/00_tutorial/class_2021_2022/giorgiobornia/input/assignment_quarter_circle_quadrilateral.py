#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/max/software/femus/applications/2021_fall/shaflynn/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("radius", 1)
notebook.set("n_q", 5)
notebook.set("half_radius", "radius*0.5")
notebook.set("n_s", 3)
notebook.set("turn_45", 45)
notebook.set("twice_radius", "2*radius")
notebook.set("back_radius", "-1*radius")
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
geomObj_1 = geompy.MakeVertex(0, 0, 0)
geomObj_2 = geompy.MakeVertex(1, 0, 0)
geomObj_3 = geompy.MakeVertex(0, 1, 0)
geomObj_4 = geompy.MakeArcCenter(geomObj_1, geomObj_3, geomObj_2,False)
geomObj_5 = geompy.MakeLineTwoPnt(geomObj_1, geomObj_3)
geomObj_6 = geompy.GetSubShape(geomObj_5, [2])
geomObj_7 = geompy.MakeLineTwoPnt(geomObj_6, geomObj_2)
geomObj_8 = geompy.MakeVertex(1, 0, 0)
geomObj_9 = geompy.MakeVertex(1, 0, 0)
geomObj_10 = geompy.MakeVertex(0.5, 0, 0)
geomObj_11 = geompy.MakeVertex(0, 1, 0)
geomObj_12 = geompy.MakeVertex(0, 0.5, 0)
geomObj_13 = geompy.MakeLineTwoPnt(O, geomObj_12)
geomObj_14 = geompy.MakeLineTwoPnt(geomObj_12, geomObj_11)
geomObj_15 = geompy.GetSubShape(geomObj_13, [2])
geomObj_16 = geompy.MakeLineTwoPnt(geomObj_10, geomObj_15)
geomObj_17 = geompy.GetSubShape(geomObj_16, [2])
geomObj_18 = geompy.MakeLineTwoPnt(geomObj_17, geomObj_9)
geomObj_19 = geompy.MakeVertex(0.5, 0.5, 0)
geomObj_20 = geompy.MakeLineTwoPnt(geomObj_12, geomObj_19)
geomObj_21 = geompy.GetSubShape(geomObj_20, [3])
geomObj_22 = geompy.MakeLineTwoPnt(geomObj_21, geomObj_17)
geomObj_23 = geompy.GetSubShape(geomObj_18, [3])
geomObj_24 = geompy.MakeArcCenter(geomObj_15, geomObj_11, geomObj_23,False)
geomObj_25 = geompy.MakeMarker(0, 0, 0, 1, 0, 0, 0, 1, 0)
geomObj_26 = geompy.MakeVertex(1, 1, 0)
geomObj_27 = geompy.GetSubShape(geomObj_14, [3])
geomObj_28 = geompy.MakeLineTwoPnt(geomObj_26, geomObj_27)
geomObj_29 = geompy.GetSubShape(geomObj_24, [3])
geomObj_30 = geompy.GetSubShape(geomObj_28, [2])
geomObj_31 = geompy.MakeLineTwoPnt(geomObj_29, geomObj_30)
geomObj_32 = geompy.MakeLineTwoPnt(geomObj_21, geomObj_30)
geomObj_33 = geompy.MakeVertexOnLinesIntersection(geomObj_32, geomObj_24)
geomObj_34 = geompy.MakeLineTwoPnt(geomObj_33, geomObj_21)
geomObj_35 = geompy.MakeFaceWires([geomObj_13, geomObj_16, geomObj_20, geomObj_22], 1)
[geomObj_36,geomObj_37,geomObj_38,geomObj_39] = geompy.ExtractShapes(geomObj_35, geompy.ShapeType["EDGE"], True)
geomObj_40 = geompy.GetSubShape(geomObj_35, [4])
[geomObj_36, geomObj_37, geomObj_38, geomObj_39, geomObj_40] = geompy.GetExistingSubObjects(geomObj_35, False)
geomObj_41 = geompy.GetSubShape(geomObj_34, [2])
geomObj_42 = geompy.MakeArcCenter(geomObj_40, geomObj_27, geomObj_41,False)
geomObj_43 = geompy.MakeArcCenter(geomObj_40, geomObj_41, geomObj_23,False)
geomObj_44 = geompy.MakeFaceWires([geomObj_14, geomObj_20, geomObj_34, geomObj_42], 1)
[geomObj_45,geomObj_46,geomObj_47,geomObj_48] = geompy.ExtractShapes(geomObj_44, geompy.ShapeType["EDGE"], True)
[geomObj_49,geomObj_50,geomObj_51,geomObj_52] = geompy.ExtractShapes(geomObj_44, geompy.ShapeType["EDGE"], True)
[geomObj_45, geomObj_46, geomObj_47, geomObj_48, geomObj_49, geomObj_50, geomObj_51, geomObj_52] = geompy.GetExistingSubObjects(geomObj_44, False)
geomObj_53 = geompy.MakeFaceWires([geomObj_18, geomObj_22, geomObj_34, geomObj_43], 1)
[geomObj_54,geomObj_55,geomObj_56,geomObj_57] = geompy.ExtractShapes(geomObj_53, geompy.ShapeType["EDGE"], True)
[geomObj_54, geomObj_55, geomObj_56, geomObj_57] = geompy.GetExistingSubObjects(geomObj_53, False)
Divided_Disk_1 = geompy.MakeDividedDisk("radius", 1, GEOM.SQUARE)
geompy.Rotate(Divided_Disk_1, OZ, "turn_45")
Cutting_Face_1 = geompy.MakeFaceHW("twice_radius", "twice_radius", 1)
Cutting_Face_2 = geompy.MakeFaceHW("twice_radius", "twice_radius", 1)
Translation_1 = geompy.MakeTranslation(Cutting_Face_1, "back_radius", 0, 0)
Translation_2 = geompy.MakeTranslation(Cutting_Face_2, 0, "back_radius", 0)
Cut_1 = geompy.MakeCutList(Divided_Disk_1, [Translation_1], True)
Cut_2 = geompy.MakeCutList(Cut_1, [Translation_2], True)
[Edge_1,Edge_2,Edge_3,Edge_4,Edge_5,Edge_6,Edge_7,Edge_8,Edge_9] = geompy.ExtractShapes(Cut_2, geompy.ShapeType["EDGE"], True)
Face_1 = geompy.MakeFaceWires([Edge_1, Edge_3, Edge_4, Edge_6], 1)
[Edge_10,Edge_11,Edge_12,Edge_13] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)
Face_2 = geompy.MakeFaceWires([Edge_2, Edge_4, Edge_5, Edge_7], 1)
[Edge_14,Edge_15,Edge_16,Edge_17] = geompy.ExtractShapes(Face_2, geompy.ShapeType["EDGE"], True)
Face_3 = geompy.MakeFaceWires([Edge_6, Edge_7, Edge_8, Edge_9], 1)
[Edge_18,Edge_19,Edge_20,Edge_21] = geompy.ExtractShapes(Face_3, geompy.ShapeType["EDGE"], True)
[Edge_10, Edge_11, Edge_12, Edge_13] = geompy.GetExistingSubObjects(Face_1, False)
[Edge_10, Edge_11, Edge_12, Edge_13] = geompy.GetExistingSubObjects(Face_1, False)
[Edge_14, Edge_15, Edge_16, Edge_17] = geompy.GetExistingSubObjects(Face_2, False)
[Edge_14, Edge_15, Edge_16, Edge_17] = geompy.GetExistingSubObjects(Face_2, False)
[Edge_18, Edge_19, Edge_20, Edge_21] = geompy.GetExistingSubObjects(Face_3, False)
[Edge_18, Edge_19, Edge_20, Edge_21] = geompy.GetExistingSubObjects(Face_3, False)
[Edge_18, Edge_19, Edge_20, Edge_21] = geompy.GetExistingSubObjects(Face_3, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Divided_Disk_1, 'Divided Disk_1' )
geompy.addToStudy( Cutting_Face_1, 'Cutting_Face_1' )
geompy.addToStudy( Cutting_Face_2, 'Cutting_Face_2' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudy( Translation_2, 'Translation_2' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudy( Cut_2, 'Cut_2' )
geompy.addToStudyInFather( Cut_2, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Cut_2, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Cut_2, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Cut_2, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( Cut_2, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Cut_2, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Cut_2, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Cut_2, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Cut_2, Edge_9, 'Edge_9' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Face_2, 'Face_2' )
geompy.addToStudy( Face_3, 'Face_3' )
geompy.addToStudyInFather( Face_1, Edge_10, 'Edge_10' )
geompy.addToStudyInFather( Face_1, Edge_11, 'Edge_11' )
geompy.addToStudyInFather( Face_1, Edge_12, 'Edge_12' )
geompy.addToStudyInFather( Face_1, Edge_13, 'Edge_13' )
geompy.addToStudyInFather( Face_2, Edge_14, 'Edge_14' )
geompy.addToStudyInFather( Face_2, Edge_15, 'Edge_15' )
geompy.addToStudyInFather( Face_2, Edge_16, 'Edge_16' )
geompy.addToStudyInFather( Face_2, Edge_17, 'Edge_17' )
geompy.addToStudyInFather( Face_3, Edge_18, 'Edge_18' )
geompy.addToStudyInFather( Face_3, Edge_19, 'Edge_19' )
geompy.addToStudyInFather( Face_3, Edge_20, 'Edge_20' )
geompy.addToStudyInFather( Face_3, Edge_21, 'Edge_21' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

smeshObj_1 = smesh.CreateHypothesisByAverageLength( 'NETGEN_Parameters_2D', 'NETGENEngine', 0.282843, 0 )
smeshObj_2 = smesh.CreateHypothesisByAverageLength( 'NETGEN_Parameters_2D', 'NETGENEngine', 0.141421, 0 )
smeshObj_3 = smesh.CreateHypothesisByAverageLength( 'NETGEN_Parameters_2D', 'NETGENEngine', 0.141421, 0 )
quad_split = smesh.CreateHypothesis('NumberOfSegments')
quad_split.SetNumberOfSegments( "n_q" )
Propagation_of_1D_Hyp = smesh.CreateHypothesis('Propagation')
Regular_1D = smesh.CreateHypothesis('Regular_1D')
Quadrangle_2D = smesh.CreateHypothesis('Quadrangle_2D')
Number_of_Segments_1 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_1.SetNumberOfSegments( 15 )
Number_of_Segments_2 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_2.SetNumberOfSegments( 15 )
Number_of_Segments_3 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_3.SetNumberOfSegments( 15 )
NETGEN_2D = smesh.CreateHypothesis('NETGEN_Remesher_2D', 'NETGENEngine')
MG_CADSurf = smesh.CreateHypothesis('MG-CADSurf_NOGEOM', 'BLSURFEngine')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
QuadFromMedialAxis_1D2D = smesh.CreateHypothesis('QuadFromMedialAxis_1D2D')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
s_Mesh_1 = smesh.Mesh(Face_1)
s_Mesh_2 = smesh.Mesh(Face_2)
Number_of_Segments_4 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_4.SetNumberOfSegments( "n_q" )
status = s_Mesh_1.AddHypothesis(Regular_1D,Edge_10)
status = s_Mesh_1.AddHypothesis(Number_of_Segments_4,Edge_10)
status = s_Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_10)
status = s_Mesh_1.AddHypothesis(Regular_1D,Edge_11)
status = s_Mesh_1.AddHypothesis(Number_of_Segments_4,Edge_11)
status = s_Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_11)
status = s_Mesh_1.AddHypothesis(Regular_1D,Edge_12)
status = s_Mesh_1.AddHypothesis(Number_of_Segments_4,Edge_12)
status = s_Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_12)
status = s_Mesh_1.AddHypothesis(Regular_1D,Edge_13)
status = s_Mesh_1.AddHypothesis(Number_of_Segments_4,Edge_13)
status = s_Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_13)
Number_of_Segments_5 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_5.SetNumberOfSegments( 15 )
status = s_Mesh_1.AddHypothesis(Regular_1D)
status = s_Mesh_1.AddHypothesis(Number_of_Segments_5)
status = s_Mesh_1.AddHypothesis(Quadrangle_2D)
isDone = s_Mesh_1.Compute()
Number_of_Segments_6 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_6.SetNumberOfSegments( "n_s" )
status = s_Mesh_2.AddHypothesis(Regular_1D,Edge_14)
status = s_Mesh_2.AddHypothesis(Regular_1D,Edge_15)
status = s_Mesh_2.AddHypothesis(Regular_1D,Edge_16)
status = s_Mesh_2.AddHypothesis(Regular_1D,Edge_17)
Number_of_Segments_7 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_7.SetNumberOfSegments( 15 )
status = s_Mesh_2.AddHypothesis(Regular_1D)
status = s_Mesh_2.AddHypothesis(Number_of_Segments_7)
status = s_Mesh_2.AddHypothesis(Quadrangle_2D)
status = s_Mesh_2.AddHypothesis(Number_of_Segments_4,Edge_14)
status = s_Mesh_2.AddHypothesis(Number_of_Segments_4,Edge_16)
status = s_Mesh_2.AddHypothesis(Number_of_Segments_4,Edge_15)
status = s_Mesh_2.AddHypothesis(Number_of_Segments_4,Edge_17)
isDone = s_Mesh_2.Compute()
Number_of_Segments_8 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_8.SetNumberOfSegments( 15 )
smesh.SetName(s_Mesh_1, 'Mesh_1')
try:
  s_Mesh_1.ExportMED(r'/home/max/software/femus/applications/2021_fall/shaflynn/input/Mesh_2_xy_assignment1_quadrangle.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
smesh.SetName(s_Mesh_1, 'Mesh_1')
try:
  s_Mesh_1.ExportMED(r'/home/max/software/femus/applications/2021_fall/shaflynn/input/Mesh_2_xy_assignment1_quadrangle.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
s_Mesh_2.ConvertToQuadratic(0, s_Mesh_2,True)
#smeshObj_4.SetName( 'Group_1_0' ) ### not created Object
#smeshObj_5.SetName( 'Group_2_0' ) ### not created Object
#smeshObj_6.SetName( 'Group_3_0' ) ### not created Object
#smeshObj_7.SetName( 'Group_4_0' ) ### not created Object
#smeshObj_8.SetName( 'Group_5_0' ) ### not created Object
#smeshObj_9.SetName( 'Group_6_0' ) ### not created Object
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
try:
  pass
except:
  #print('ExportMED() failed. Invalid file name?') ### not created Object
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
s_Mesh_3 = smesh.Mesh(Face_3)
Edge_18_1 = s_Mesh_3.GroupOnGeom(Edge_18,'Edge_18',SMESH.EDGE)
Edge_19_1 = s_Mesh_3.GroupOnGeom(Edge_19,'Edge_19',SMESH.EDGE)
Edge_20_1 = s_Mesh_3.GroupOnGeom(Edge_20,'Edge_20',SMESH.EDGE)
Edge_21_1 = s_Mesh_3.GroupOnGeom(Edge_21,'Edge_21',SMESH.EDGE)
status = s_Mesh_3.AddHypothesis(Regular_1D,Edge_18)
status = s_Mesh_3.AddHypothesis(Number_of_Segments_4,Edge_18)
status = s_Mesh_3.AddHypothesis(Regular_1D,Edge_19)
status = s_Mesh_3.AddHypothesis(Number_of_Segments_4,Edge_19)
status = s_Mesh_3.AddHypothesis(Regular_1D,Edge_20)
status = s_Mesh_3.AddHypothesis(Number_of_Segments_4,Edge_20)
status = s_Mesh_3.AddHypothesis(Regular_1D,Edge_21)
status = s_Mesh_3.AddHypothesis(Number_of_Segments_4,Edge_21)
Number_of_Segments_9 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_9.SetNumberOfSegments( 15 )
status = s_Mesh_3.AddHypothesis(Regular_1D)
status = s_Mesh_3.AddHypothesis(Number_of_Segments_9)
status = s_Mesh_3.AddHypothesis(Quadrangle_2D)
isDone = s_Mesh_3.Compute()
[ Edge_18_1, Edge_19_1, Edge_20_1, Edge_21_1 ] = s_Mesh_3.GetGroups()
s_Mesh_3.ConvertToQuadratic(0, s_Mesh_3,True)
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Mesh_1 = smesh.Concatenate( [ s_Mesh_1.GetMesh(), s_Mesh_2.GetMesh(), s_Mesh_3.GetMesh() ], 1, 1, 1e-05, False )
[ Edge_18_2, Edge_19_2, Edge_20_2, Edge_21_2 ] = Mesh_1.GetGroups()
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.EDGE,SMESH.FT_FreeBorders,SMESH.FT_Undefined,0,SMESH.FT_LogicalNOT)
aCriteria.append(aCriterion)
#isDone = smeshObj_10.RemoveElements( [ 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 51, 52, 53, 54, 55 ] ) ### not created Object
[ Edge_18_2, Edge_19_2, Edge_20_2, Edge_21_2 ] = Mesh_1.GetGroups()
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_1_0' )
nbAdd = Group_1_0.Add( [ 1, 2, 3, 4, 5, 36, 37, 38, 39, 40 ] )
[ Edge_18_2, Edge_19_2, Edge_20_2, Edge_21_2, Group_1_0 ] = Mesh_1.GetGroups()
Group_2_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_2_0' )
nbAdd = Group_2_0.Add( [ 41, 42, 43, 44, 45, 71, 72, 73, 74, 75 ] )
[ Edge_18_2, Edge_19_2, Edge_20_2, Edge_21_2, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Group_3_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_3_0' )
nbAdd = Group_3_0.Add( [ 6, 7, 8, 9, 10, 76, 77, 78, 79, 80 ] )
[ Edge_18_2, Edge_19_2, Edge_20_2, Edge_21_2, Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
isDone = Mesh_1.Compute()
[ Edge_18_2, Edge_19_2, Edge_20_2, Edge_21_2, Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
Sub_mesh_1 = s_Mesh_1.GetSubMesh( Edge_10, 'Sub-mesh_1' )
Sub_mesh_2 = s_Mesh_1.GetSubMesh( Edge_11, 'Sub-mesh_2' )
Sub_mesh_3 = s_Mesh_1.GetSubMesh( Edge_12, 'Sub-mesh_3' )
Sub_mesh_4 = s_Mesh_1.GetSubMesh( Edge_13, 'Sub-mesh_4' )
Sub_mesh_5 = s_Mesh_2.GetSubMesh( Edge_14, 'Sub-mesh_5' )
Sub_mesh_6 = s_Mesh_2.GetSubMesh( Edge_15, 'Sub-mesh_6' )
Sub_mesh_7 = s_Mesh_2.GetSubMesh( Edge_16, 'Sub-mesh_7' )
Sub_mesh_8 = s_Mesh_2.GetSubMesh( Edge_17, 'Sub-mesh_8' )
Sub_mesh_9 = s_Mesh_3.GetSubMesh( Edge_18, 'Sub-mesh_9' )
Sub_mesh_10 = s_Mesh_3.GetSubMesh( Edge_19, 'Sub-mesh_10' )
Sub_mesh_11 = s_Mesh_3.GetSubMesh( Edge_20, 'Sub-mesh_11' )
Sub_mesh_12 = s_Mesh_3.GetSubMesh( Edge_21, 'Sub-mesh_12' )

## some objects were removed
aStudyBuilder = salome.myStudy.NewBuilder()
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_3))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_2))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_1))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)

## Set names of Mesh objects
smesh.SetName(Number_of_Segments_7, 'Number of Segments_7')
smesh.SetName(Number_of_Segments_8, 'Number of Segments_8')
smesh.SetName(Number_of_Segments_5, 'Number of Segments_5')
smesh.SetName(Number_of_Segments_6, 'Number of Segments_6')
smesh.SetName(Number_of_Segments_9, 'Number of Segments_9')
smesh.SetName(Sub_mesh_8, 'Sub-mesh_8')
smesh.SetName(s_Mesh_3.GetMesh(), 's_Mesh_3')
smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Regular_1D, 'Regular_1D')
smesh.SetName(Sub_mesh_7, 'Sub-mesh_7')
smesh.SetName(Quadrangle_2D, 'Quadrangle_2D')
smesh.SetName(Sub_mesh_6, 'Sub-mesh_6')
smesh.SetName(NETGEN_2D, 'NETGEN 2D')
smesh.SetName(MG_CADSurf, 'MG-CADSurf')
smesh.SetName(QuadFromMedialAxis_1D2D, 'QuadFromMedialAxis_1D2D')
smesh.SetName(s_Mesh_2.GetMesh(), 's_Mesh_2')
smesh.SetName(s_Mesh_1.GetMesh(), 's_Mesh_1')
smesh.SetName(Sub_mesh_9, 'Sub-mesh_9')
smesh.SetName(Edge_21_2, 'Edge_21')
smesh.SetName(Sub_mesh_10, 'Sub-mesh_10')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Sub_mesh_11, 'Sub-mesh_11')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Sub_mesh_12, 'Sub-mesh_12')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Edge_18_2, 'Edge_18')
smesh.SetName(Edge_19_2, 'Edge_19')
smesh.SetName(Edge_20_2, 'Edge_20')
smesh.SetName(Edge_21_1, 'Edge_21')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Edge_19_1, 'Edge_19')
smesh.SetName(Edge_20_1, 'Edge_20')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Edge_18_1, 'Edge_18')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(quad_split, 'quad_split')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
