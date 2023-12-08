#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.8.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/applications/2021_fall/hdongamm/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("n_x", 3)
notebook.set("n_y", 6)
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
Divided_Disk_1 = geompy.MakeDividedDisk(1, 1, GEOM.SQUARE)
geompy.Rotate(Divided_Disk_1, OZ, 45*math.pi/180.0)
Face_1 = geompy.MakeFaceHW(2, 2, 1)
Translation_1 = geompy.MakeTranslation(Face_1, 0, -1, 0)
Cut_1 = geompy.MakeCutList(Divided_Disk_1, [Translation_1], True)
[Face_2,Face_3,Face_4,Face_5] = geompy.ExtractShapes(Cut_1, geompy.ShapeType["FACE"], True)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_2, geompy.ShapeType["EDGE"], True)
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Face_3, geompy.ShapeType["EDGE"], True)
[Edge_9,Edge_10,Edge_11,Edge_12] = geompy.ExtractShapes(Face_4, geompy.ShapeType["EDGE"], True)
[Edge_13,Edge_14,Edge_15,Edge_16] = geompy.ExtractShapes(Face_5, geompy.ShapeType["EDGE"], True)
[Face_2, Face_3, Face_4, Face_5] = geompy.GetExistingSubObjects(Cut_1, False)
Auto_group_for_Group_1_0 = geompy.CreateGroup(Cut_1, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_1_0, [Edge_2, Edge_6, Edge_15])
Auto_group_for_Group_2_0 = geompy.CreateGroup(Cut_1, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_2_0, [Edge_1, Edge_11, Edge_16])
O_1 = geompy.MakeVertex(0, 0, 0)
OX_1 = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY_1 = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ_1 = geompy.MakeVectorDXDYDZ(0, 0, 1)
Divided_Disk_1_1 = geompy.MakeDividedDisk(1, 1, GEOM.SQUARE)
geompy.Rotate(Divided_Disk_1_1, OZ_1, 45*math.pi/180.0)
Face_1_1 = geompy.MakeFaceHW(2, 2, 1)
Translation_1_1 = geompy.MakeTranslation(Face_1_1, 0, -1, 0)
Cut_1_1 = geompy.MakeCutList(Divided_Disk_1_1, [Translation_1_1], True)
[Face_2_1,Face_3_1,Face_4_1,Face_5_1] = geompy.ExtractShapes(Cut_1_1, geompy.ShapeType["FACE"], True)
[Edge_1_1,Edge_2_1,Edge_3_1,Edge_4_1] = geompy.ExtractShapes(Face_2_1, geompy.ShapeType["EDGE"], True)
[Edge_5_1,Edge_6_1,Edge_7_1,Edge_8_1] = geompy.ExtractShapes(Face_3_1, geompy.ShapeType["EDGE"], True)
[Edge_9_1,Edge_10_1,Edge_11_1,Edge_12_1] = geompy.ExtractShapes(Face_4_1, geompy.ShapeType["EDGE"], True)
[Edge_13_1,Edge_14_1,Edge_15_1,Edge_16_1] = geompy.ExtractShapes(Face_5_1, geompy.ShapeType["EDGE"], True)
[Face_2_1, Face_3_1, Face_4_1, Face_5_1] = geompy.GetExistingSubObjects(Cut_1_1, False)
Auto_group_for_Group_1_0_1 = geompy.CreateGroup(Cut_1_1, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_1_0_1, [Edge_2_1, Edge_6_1, Edge_15_1])
Auto_group_for_Group_2_0_1 = geompy.CreateGroup(Cut_1_1, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_2_0_1, [Edge_1_1, Edge_11_1, Edge_16_1])
O_2 = geompy.MakeVertex(0, 0, 0)
OX_2 = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY_2 = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ_2 = geompy.MakeVectorDXDYDZ(0, 0, 1)
Divided_Disk_1_2 = geompy.MakeDividedDisk(1, 1, GEOM.SQUARE)
geompy.Rotate(Divided_Disk_1_2, OZ_2, 45*math.pi/180.0)
Face_1_2 = geompy.MakeFaceHW(2, 2, 1)
Translation_1_2 = geompy.MakeTranslation(Face_1_2, 0, -1, 0)
Cut_1_2 = geompy.MakeCutList(Divided_Disk_1_2, [Translation_1_2], True)
[Face_2_2,Face_3_2,Face_4_2,Face_5_2] = geompy.ExtractShapes(Cut_1_2, geompy.ShapeType["FACE"], True)
[Edge_1_2,Edge_2_2,Edge_3_2,Edge_4_2] = geompy.ExtractShapes(Face_2_2, geompy.ShapeType["EDGE"], True)
[Edge_5_2,Edge_6_2,Edge_7_2,Edge_8_2] = geompy.ExtractShapes(Face_3_2, geompy.ShapeType["EDGE"], True)
[Edge_9_2,Edge_10_2,Edge_11_2,Edge_12_2] = geompy.ExtractShapes(Face_4_2, geompy.ShapeType["EDGE"], True)
[Edge_13_2,Edge_14_2,Edge_15_2,Edge_16_2] = geompy.ExtractShapes(Face_5_2, geompy.ShapeType["EDGE"], True)
[Face_2_2, Face_3_2, Face_4_2, Face_5_2] = geompy.GetExistingSubObjects(Cut_1_2, False)
Auto_group_for_Group_1_0_2 = geompy.CreateGroup(Cut_1_2, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_1_0_2, [Edge_2_2, Edge_6_2, Edge_15_2])
Auto_group_for_Group_2_0_2 = geompy.CreateGroup(Cut_1_2, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_2_0_2, [Edge_1_2, Edge_11_2, Edge_16_2])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Divided_Disk_1, 'Divided Disk_1' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudyInFather( Cut_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Cut_1, Face_3, 'Face_3' )
geompy.addToStudyInFather( Cut_1, Face_4, 'Face_4' )
geompy.addToStudyInFather( Cut_1, Face_5, 'Face_5' )
geompy.addToStudyInFather( Face_2, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_2, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_2, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Face_2, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( Face_3, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_3, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_3, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Face_3, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Face_4, Edge_9, 'Edge_9' )
geompy.addToStudyInFather( Face_4, Edge_10, 'Edge_10' )
geompy.addToStudyInFather( Face_4, Edge_11, 'Edge_11' )
geompy.addToStudyInFather( Face_4, Edge_12, 'Edge_12' )
geompy.addToStudyInFather( Face_5, Edge_13, 'Edge_13' )
geompy.addToStudyInFather( Face_5, Edge_14, 'Edge_14' )
geompy.addToStudyInFather( Face_5, Edge_15, 'Edge_15' )
geompy.addToStudyInFather( Face_5, Edge_16, 'Edge_16' )
geompy.addToStudyInFather( Cut_1, Auto_group_for_Group_1_0, 'Auto_group_for_Group_1_0' )
geompy.addToStudyInFather( Cut_1, Auto_group_for_Group_2_0, 'Auto_group_for_Group_2_0' )
geompy.addToStudy( O_1, 'O' )
geompy.addToStudy( OX_1, 'OX' )
geompy.addToStudy( OY_1, 'OY' )
geompy.addToStudy( OZ_1, 'OZ' )
geompy.addToStudy( Divided_Disk_1_1, 'Divided Disk_1' )
geompy.addToStudy( Face_1_1, 'Face_1' )
geompy.addToStudy( Translation_1_1, 'Translation_1' )
geompy.addToStudy( Cut_1_1, 'Cut_1' )
geompy.addToStudyInFather( Cut_1_1, Face_2_1, 'Face_2' )
geompy.addToStudyInFather( Cut_1_1, Face_3_1, 'Face_3' )
geompy.addToStudyInFather( Cut_1_1, Face_4_1, 'Face_4' )
geompy.addToStudyInFather( Cut_1_1, Face_5_1, 'Face_5' )
geompy.addToStudyInFather( Face_2_1, Edge_1_1, 'Edge_1' )
geompy.addToStudyInFather( Face_2_1, Edge_2_1, 'Edge_2' )
geompy.addToStudyInFather( Face_2_1, Edge_3_1, 'Edge_3' )
geompy.addToStudyInFather( Face_2_1, Edge_4_1, 'Edge_4' )
geompy.addToStudyInFather( Face_3_1, Edge_5_1, 'Edge_5' )
geompy.addToStudyInFather( Face_3_1, Edge_6_1, 'Edge_6' )
geompy.addToStudyInFather( Face_3_1, Edge_7_1, 'Edge_7' )
geompy.addToStudyInFather( Face_3_1, Edge_8_1, 'Edge_8' )
geompy.addToStudyInFather( Face_4_1, Edge_9_1, 'Edge_9' )
geompy.addToStudyInFather( Face_4_1, Edge_10_1, 'Edge_10' )
geompy.addToStudyInFather( Face_4_1, Edge_11_1, 'Edge_11' )
geompy.addToStudyInFather( Face_4_1, Edge_12_1, 'Edge_12' )
geompy.addToStudyInFather( Face_5_1, Edge_13_1, 'Edge_13' )
geompy.addToStudyInFather( Face_5_1, Edge_14_1, 'Edge_14' )
geompy.addToStudyInFather( Face_5_1, Edge_15_1, 'Edge_15' )
geompy.addToStudyInFather( Face_5_1, Edge_16_1, 'Edge_16' )
geompy.addToStudyInFather( Cut_1_1, Auto_group_for_Group_1_0_1, 'Auto_group_for_Group_1_0' )
geompy.addToStudyInFather( Cut_1_1, Auto_group_for_Group_2_0_1, 'Auto_group_for_Group_2_0' )
geompy.addToStudy( O_2, 'O' )
geompy.addToStudy( OX_2, 'OX' )
geompy.addToStudy( OY_2, 'OY' )
geompy.addToStudy( OZ_2, 'OZ' )
geompy.addToStudy( Divided_Disk_1_2, 'Divided Disk_1' )
geompy.addToStudy( Face_1_2, 'Face_1' )
geompy.addToStudy( Translation_1_2, 'Translation_1' )
geompy.addToStudy( Cut_1_2, 'Cut_1' )
geompy.addToStudyInFather( Cut_1_2, Face_2_2, 'Face_2' )
geompy.addToStudyInFather( Cut_1_2, Face_3_2, 'Face_3' )
geompy.addToStudyInFather( Cut_1_2, Face_4_2, 'Face_4' )
geompy.addToStudyInFather( Cut_1_2, Face_5_2, 'Face_5' )
geompy.addToStudyInFather( Face_2_2, Edge_1_2, 'Edge_1' )
geompy.addToStudyInFather( Face_2_2, Edge_2_2, 'Edge_2' )
geompy.addToStudyInFather( Face_2_2, Edge_3_2, 'Edge_3' )
geompy.addToStudyInFather( Face_2_2, Edge_4_2, 'Edge_4' )
geompy.addToStudyInFather( Face_3_2, Edge_5_2, 'Edge_5' )
geompy.addToStudyInFather( Face_3_2, Edge_6_2, 'Edge_6' )
geompy.addToStudyInFather( Face_3_2, Edge_7_2, 'Edge_7' )
geompy.addToStudyInFather( Face_3_2, Edge_8_2, 'Edge_8' )
geompy.addToStudyInFather( Face_4_2, Edge_9_2, 'Edge_9' )
geompy.addToStudyInFather( Face_4_2, Edge_10_2, 'Edge_10' )
geompy.addToStudyInFather( Face_4_2, Edge_11_2, 'Edge_11' )
geompy.addToStudyInFather( Face_4_2, Edge_12_2, 'Edge_12' )
geompy.addToStudyInFather( Face_5_2, Edge_13_2, 'Edge_13' )
geompy.addToStudyInFather( Face_5_2, Edge_14_2, 'Edge_14' )
geompy.addToStudyInFather( Face_5_2, Edge_15_2, 'Edge_15' )
geompy.addToStudyInFather( Face_5_2, Edge_16_2, 'Edge_16' )
geompy.addToStudyInFather( Cut_1_2, Auto_group_for_Group_1_0_2, 'Auto_group_for_Group_1_0' )
geompy.addToStudyInFather( Cut_1_2, Auto_group_for_Group_2_0_2, 'Auto_group_for_Group_2_0' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Cut_1)
Regular_1D = Mesh_1.Segment(geom=Edge_1)
Number_of_Segments_1 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_1.SetNumberOfSegments( "n_x" )
Number_of_Segments_1.SetReversedEdges( [] )
Number_of_Segments_1.SetObjectEntry( "Cut_1" )
Regular_1D_1 = Mesh_1.Segment(geom=Edge_2)
Number_of_Segments_2 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_2.SetNumberOfSegments( "n_x" )
Number_of_Segments_2.SetReversedEdges( [] )
Number_of_Segments_2.SetObjectEntry( "Cut_1" )
Regular_1D_2 = Mesh_1.Segment(geom=Edge_6)
Number_of_Segments_3 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_3.SetNumberOfSegments( "n_y" )
Number_of_Segments_3.SetReversedEdges( [] )
Number_of_Segments_3.SetObjectEntry( "Cut_1" )
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Number_of_Segments_4 = Regular_1D.NumberOfSegments("n_x",None,[])
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Number_of_Segments_5 = Regular_1D_1.NumberOfSegments("n_x",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_2)
Number_of_Segments_6 = Regular_1D_2.NumberOfSegments("n_y",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_6)
Regular_1D_3 = Mesh_1.Segment()
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Group_1_0 = Mesh_1.GroupOnGeom(Auto_group_for_Group_1_0,'Group_1_0',SMESH.EDGE)
Group_2_0 = Mesh_1.GroupOnGeom(Auto_group_for_Group_2_0,'Group_2_0',SMESH.EDGE)
isDone = Mesh_1.RemoveElements( [ 1, 2, 3, 10, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 28, 29, 30 ] )
Regular_1D_4 = smesh.CreateHypothesis( "Regular_1D" )
NumberOfSegments_n_x = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_x.SetNumberOfSegments( "n_x" )
NumberOfSegments_n_x.SetReversedEdges( [] )
NumberOfSegments_n_x.SetObjectEntry( "Cut_1_1" )
Regular_1D_5 = smesh.CreateHypothesis( "Regular_1D" )
NumberOfSegments_n_x_1 = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_x_1.SetNumberOfSegments( "n_x" )
NumberOfSegments_n_x_1.SetReversedEdges( [] )
NumberOfSegments_n_x_1.SetObjectEntry( "Cut_1_1" )
Regular_1D_6 = smesh.CreateHypothesis( "Regular_1D" )
NumberOfSegments_n_y = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_y.SetNumberOfSegments( "n_y" )
NumberOfSegments_n_y.SetReversedEdges( [] )
NumberOfSegments_n_y.SetObjectEntry( "Cut_1_1" )
Quadrangle_2D_1 = smesh.CreateHypothesis( "Quadrangle_2D" )
NumberOfSegments_n_x_2 = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_x_2.SetNumberOfSegments( "n_x" )
NumberOfSegments_n_x_2.SetReversedEdges( [] )
NumberOfSegments_n_x_2.SetObjectEntry( "Cut_1_1" )
NumberOfSegments_n_x_3 = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_x_3.SetNumberOfSegments( "n_x" )
NumberOfSegments_n_x_3.SetReversedEdges( [] )
NumberOfSegments_n_x_3.SetObjectEntry( "Cut_1_1" )
NumberOfSegments_n_y_1 = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_y_1.SetNumberOfSegments( "n_y" )
NumberOfSegments_n_y_1.SetReversedEdges( [] )
NumberOfSegments_n_y_1.SetObjectEntry( "Cut_1_1" )
Regular_1D_7 = smesh.CreateHypothesis( "Regular_1D" )
Regular_1D_8 = smesh.CreateHypothesis( "Regular_1D" )
NumberOfSegments_n_x_4 = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_x_4.SetNumberOfSegments( "n_x" )
NumberOfSegments_n_x_4.SetReversedEdges( [] )
NumberOfSegments_n_x_4.SetObjectEntry( "Cut_1_2" )
Regular_1D_9 = smesh.CreateHypothesis( "Regular_1D" )
NumberOfSegments_n_x_5 = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_x_5.SetNumberOfSegments( "n_x" )
NumberOfSegments_n_x_5.SetReversedEdges( [] )
NumberOfSegments_n_x_5.SetObjectEntry( "Cut_1_2" )
Regular_1D_10 = smesh.CreateHypothesis( "Regular_1D" )
NumberOfSegments_n_y_2 = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_y_2.SetNumberOfSegments( "n_y" )
NumberOfSegments_n_y_2.SetReversedEdges( [] )
NumberOfSegments_n_y_2.SetObjectEntry( "Cut_1_2" )
Quadrangle_2D_2 = smesh.CreateHypothesis( "Quadrangle_2D" )
NumberOfSegments_n_x_6 = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_x_6.SetNumberOfSegments( "n_x" )
NumberOfSegments_n_x_6.SetReversedEdges( [] )
NumberOfSegments_n_x_6.SetObjectEntry( "Cut_1_2" )
NumberOfSegments_n_x_7 = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_x_7.SetNumberOfSegments( "n_x" )
NumberOfSegments_n_x_7.SetReversedEdges( [] )
NumberOfSegments_n_x_7.SetObjectEntry( "Cut_1_2" )
NumberOfSegments_n_y_3 = smesh.CreateHypothesis('NumberOfSegments')
NumberOfSegments_n_y_3.SetNumberOfSegments( "n_y" )
NumberOfSegments_n_y_3.SetReversedEdges( [] )
NumberOfSegments_n_y_3.SetObjectEntry( "Cut_1_2" )
Regular_1D_11 = smesh.CreateHypothesis( "Regular_1D" )
Group_2_0.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))
Group_1_0.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))
Sub_mesh_1 = Regular_1D.GetSubMesh()
Sub_mesh_2 = Regular_1D_1.GetSubMesh()
Sub_mesh_3 = Regular_1D_2.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_6')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_5, 'Number of Segments_5')
smesh.SetName(Number_of_Segments_6, 'Number of Segments_6')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(NumberOfSegments_n_x, 'NumberOfSegments=n_x,[],0:1:1:16')
smesh.SetName(NumberOfSegments_n_x_1, 'NumberOfSegments=n_x,[],0:1:1:16')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(NumberOfSegments_n_y_1, 'NumberOfSegments=n_y,[],0:1:1:16')
smesh.SetName(NumberOfSegments_n_x_3, 'NumberOfSegments=n_x,[],0:1:1:16')
smesh.SetName(NumberOfSegments_n_x_2, 'NumberOfSegments=n_x,[],0:1:1:16')
smesh.SetName(NumberOfSegments_n_y, 'NumberOfSegments=n_y,[],0:1:1:16')
smesh.SetName(NumberOfSegments_n_x_6, 'NumberOfSegments=n_x,[],0:1:1:24')
smesh.SetName(NumberOfSegments_n_y_2, 'NumberOfSegments=n_y,[],0:1:1:24')
smesh.SetName(NumberOfSegments_n_x_5, 'NumberOfSegments=n_x,[],0:1:1:24')
smesh.SetName(NumberOfSegments_n_x_4, 'NumberOfSegments=n_x,[],0:1:1:24')
smesh.SetName(NumberOfSegments_n_y_3, 'NumberOfSegments=n_y,[],0:1:1:24')
smesh.SetName(NumberOfSegments_n_x_7, 'NumberOfSegments=n_x,[],0:1:1:24')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()