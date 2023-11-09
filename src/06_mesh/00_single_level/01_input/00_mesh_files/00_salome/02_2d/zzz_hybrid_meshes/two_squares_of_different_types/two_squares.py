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
Face_1 = geompy.MakeFaceHW(1, 2, 1)
geompy.TranslateDXDYDZ(Face_1, 0.5, 0, 0)
Partition_1 = geompy.MakePartition([Face_1], [OX], [], [], geompy.ShapeType["FACE"], 0, [], 1)
[geomObj_1,geomObj_2,geomObj_3,geomObj_4,geomObj_5,geomObj_6,geomObj_7] = geompy.ExtractShapes(Partition_1, geompy.ShapeType["EDGE"], True)
[Face_2,Face_3] = geompy.ExtractShapes(Partition_1, geompy.ShapeType["FACE"], True)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_2, geompy.ShapeType["EDGE"], True)
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Face_3, geompy.ShapeType["EDGE"], True)
Auto_group_for_Group_1_0 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_1_0, [Edge_1, Edge_2, Edge_5])
Auto_group_for_Group_2 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Group_2, [Edge_4, Edge_7, Edge_8])
Auto_group_for_Group_1 = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
geompy.UnionList(Auto_group_for_Group_1, [Face_2, Face_3])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Partition_1, 'Partition_1' )
geompy.addToStudyInFather( Partition_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Partition_1, Face_3, 'Face_3' )
geompy.addToStudyInFather( Face_2, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_2, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_2, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Face_2, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( Face_3, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_3, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_3, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Face_3, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Partition_1, Auto_group_for_Group_1_0, 'Auto_group_for_Group_1_0' )
geompy.addToStudyInFather( Partition_1, Auto_group_for_Group_2, 'Auto_group_for_Group_2' )
geompy.addToStudyInFather( Partition_1, Auto_group_for_Group_1, 'Auto_group_for_Group_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Number_of_Segments_1 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_1.SetNumberOfSegments( 1 )
Number_of_Segments_2 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_2.SetNumberOfSegments( 1 )
Propagation_of_1D_Hyp = smesh.CreateHypothesis('Propagation')
Mesh_1 = smesh.Mesh(Partition_1,'Mesh_1')
Regular_1D = Mesh_1.Segment(geom=Face_2)
Sub_mesh_1 = Regular_1D.GetSubMesh()
Number_of_Segments_3 = Regular_1D.NumberOfSegments(1)
Propagation_of_1D_Hyp_1 = Regular_1D.Propagation()
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE,geom=Face_2)
Regular_1D_1 = Mesh_1.Segment(geom=Face_3)
Sub_mesh_2 = Regular_1D_1.GetSubMesh()
Number_of_Segments_4 = Regular_1D_1.NumberOfSegments(1)
Propagation_of_1D_Hyp_2 = Regular_1D_1.Propagation()
Quadrangle_2D_1 = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE,geom=Face_3)
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1, Sub_mesh_2 ] ])
Group_1_0 = Mesh_1.GroupOnGeom(Auto_group_for_Group_1_0,'Group_1_0',SMESH.EDGE)
Group_2_0 = Mesh_1.GroupOnGeom(Auto_group_for_Group_2,'Group_2',SMESH.EDGE)
Group_2_0.SetName( 'Group_2_0' )
Group_6_0 = Mesh_1.GroupOnGeom(Auto_group_for_Group_1,'Group_1',SMESH.FACE)
isDone = Mesh_1.Compute()
[ Group_1_0, Group_2_0, Group_6_0 ] = Mesh_1.GetGroups()
Group_6_0.SetName( 'Group_6_0' )
[ Group_1_0, Group_2_0, Group_6_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0, Sub_mesh_1)
Mesh_1.ConvertToQuadratic(0, Sub_mesh_2,True)
Group_6_0.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))
Group_2_0.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))
Group_1_0.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Propagation_of_1D_Hyp_2, 'Propagation of 1D Hyp. on Opposite Edges_3')
smesh.SetName(Group_6_0, 'Group_6_0')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Propagation_of_1D_Hyp_1, 'Propagation of 1D Hyp. on Opposite Edges_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
