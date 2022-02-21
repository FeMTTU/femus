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
[Edge_14, Edge_15, Edge_16, Edge_17] = geompy.GetExistingSubObjects(Face_2, False)
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
Mesh_1 = smesh.Mesh(Face_1)
Regular_1D = Mesh_1.Segment(geom=Edge_10)
quad_split = Regular_1D.NumberOfSegments("n_q")
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Regular_1D_1 = Mesh_1.Segment(geom=Edge_11)
status = Mesh_1.AddHypothesis(quad_split,Edge_11)
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_11)
Regular_1D_2 = Mesh_1.Segment(geom=Edge_12)
status = Mesh_1.AddHypothesis(quad_split,Edge_12)
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_12)
Regular_1D_3 = Mesh_1.Segment(geom=Edge_13)
status = Mesh_1.AddHypothesis(quad_split,Edge_13)
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_13)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_4 = Mesh_1.Segment()
Number_of_Segments_1 = Regular_1D_4.NumberOfSegments(15)
isDone = Mesh_1.Compute()
Regular_1D_5 = smesh.CreateHypothesis( "Regular_1D" )
Regular_1D_6 = smesh.CreateHypothesis( "Regular_1D" )
Regular_1D_7 = smesh.CreateHypothesis( "Regular_1D" )
Regular_1D_8 = smesh.CreateHypothesis( "Regular_1D" )
Number_of_Segments_2 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_2.SetNumberOfSegments( 15 )
Regular_1D_9 = smesh.CreateHypothesis( "Regular_1D" )
Quadrangle_2D_1 = smesh.CreateHypothesis( "Quadrangle_2D" )
Number_of_Segments_3 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_3.SetNumberOfSegments( 15 )
Regular_1D_10 = smesh.CreateHypothesis( "Regular_1D" )
Quadrangle_2D_2 = smesh.CreateHypothesis( "Quadrangle_2D" )
Regular_1D_11 = smesh.CreateHypothesis( "Regular_1D" )
Regular_1D_12 = smesh.CreateHypothesis( "Regular_1D" )
Regular_1D_13 = smesh.CreateHypothesis( "Regular_1D" )
Regular_1D_14 = smesh.CreateHypothesis( "Regular_1D" )
Mesh_1.ConvertToQuadratic(0)
Mesh_1.ConvertFromQuadratic()
Mesh_1.ConvertToQuadratic(0)
NETGEN_2D = smesh.CreateHypothesis('NETGEN_Remesher_2D', 'NETGENEngine')
MG_CADSurf = smesh.CreateHypothesis('MG-CADSurf_NOGEOM', 'BLSURFEngine')
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/max/software/femus/applications/2021_fall/shaflynn/input/Mesh_1_assignment_1_quadrangle.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
isDone = Mesh_1.Compute()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/max/software/femus/applications/2021_fall/shaflynn/input/Mesh_1_assignment_1_quadrangle.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/max/software/femus/applications/2021_fall/shaflynn/input/Mesh_1_assignment_1_quadrangle.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
status = Mesh_1.RemoveHypothesis(Quadrangle_2D)
QuadFromMedialAxis_1D2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUAD_MA_PROJ)
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0)
Group_1_0 = Mesh_1.GroupOnGeom(Edge_10,'Group_1_0',SMESH.EDGE)
[ Group_1_0 ] = Mesh_1.GetGroups()
Group_2_0 = Mesh_1.GroupOnGeom(Edge_11,'Group_2_0',SMESH.EDGE)
[ Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Group_3_0 = Mesh_1.GroupOnGeom(Edge_12,'Group_3_0',SMESH.EDGE)
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
Group_4_0 = Mesh_1.GroupOnGeom(Edge_13,'Group_4_0',SMESH.EDGE)
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0 ] = Mesh_1.GetGroups()
isDone = Mesh_1.Compute()
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/max/software/femus/applications/2021_fall/shaflynn/input/Mesh_1_assignment_1_quadrangle.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Sub_mesh_1 = Regular_1D.GetSubMesh()
Sub_mesh_2 = Regular_1D_1.GetSubMesh()
Sub_mesh_3 = Regular_1D_2.GetSubMesh()
Sub_mesh_4 = Regular_1D_3.GetSubMesh()

## some objects were removed
aStudyBuilder = salome.myStudy.NewBuilder()
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_3))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_2))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_1))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)

## Set names of Mesh objects
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(NETGEN_2D, 'NETGEN 2D')
smesh.SetName(MG_CADSurf, 'MG-CADSurf')
smesh.SetName(QuadFromMedialAxis_1D2D.GetAlgorithm(), 'QuadFromMedialAxis_1D2D')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(quad_split, 'quad_split')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
