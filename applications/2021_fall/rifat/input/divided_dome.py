#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.3.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/applications/OptimalControl/00_laplacian_all_dims/input')

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
Sphere_1 = geompy.MakeSphereR(1)
[Face_1] = geompy.ExtractShapes(Sphere_1, geompy.ShapeType["FACE"], True)
Divided_Cylinder_1 = geompy.MakeDividedCylinder(1, 2, GEOM.SQUARE)
Common_1 = geompy.MakeCommonList([Divided_Cylinder_1, Face_1], True)
[Face_2,Face_3,Face_4,Face_5,Face_6] = geompy.ExtractShapes(Common_1, geompy.ShapeType["FACE"], True)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_2, geompy.ShapeType["EDGE"], True)
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Face_3, geompy.ShapeType["EDGE"], True)
[Edge_9,Edge_10,Edge_11,Edge_12] = geompy.ExtractShapes(Face_5, geompy.ShapeType["EDGE"], True)
[Edge_13,Edge_14,Edge_15,Edge_16] = geompy.ExtractShapes(Face_6, geompy.ShapeType["EDGE"], True)
Thickness_1 = geompy.MakeThickSolid(Face_4, 1, [])
NoExtraEdges_1 = geompy.RemoveExtraEdges(Thickness_1, False)
[Face_7,Face_8,Face_9,Face_10,Face_11,Face_12] = geompy.ExtractShapes(NoExtraEdges_1, geompy.ShapeType["FACE"], True)
[Edge_17,Edge_18,Edge_19,Edge_20] = geompy.ExtractShapes(Face_9, geompy.ShapeType["EDGE"], True)
Shell_1 = geompy.MakeShell([Face_2, Face_3, Face_5, Face_6, Face_9])
[Edge_21,Edge_22,Edge_23,Edge_24,Edge_25,Edge_26,Edge_27,Edge_28,Edge_29,Edge_30,Edge_31,Edge_32] = geompy.ExtractShapes(Shell_1, geompy.ShapeType["EDGE"], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Sphere_1, 'Sphere_1' )
geompy.addToStudyInFather( Sphere_1, Face_1, 'Face_1' )
geompy.addToStudy( Divided_Cylinder_1, 'Divided Cylinder_1' )
geompy.addToStudy( Common_1, 'Common_1' )
geompy.addToStudyInFather( Common_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Common_1, Face_3, 'Face_3' )
geompy.addToStudyInFather( Common_1, Face_4, 'Face_4' )
geompy.addToStudyInFather( Common_1, Face_5, 'Face_5' )
geompy.addToStudyInFather( Common_1, Face_6, 'Face_6' )
geompy.addToStudyInFather( Face_2, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_2, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_2, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Face_2, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( Face_3, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_3, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_3, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Face_3, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Face_5, Edge_9, 'Edge_9' )
geompy.addToStudyInFather( Face_5, Edge_10, 'Edge_10' )
geompy.addToStudyInFather( Face_5, Edge_11, 'Edge_11' )
geompy.addToStudyInFather( Face_5, Edge_12, 'Edge_12' )
geompy.addToStudyInFather( Face_6, Edge_13, 'Edge_13' )
geompy.addToStudyInFather( Face_6, Edge_14, 'Edge_14' )
geompy.addToStudyInFather( Face_6, Edge_15, 'Edge_15' )
geompy.addToStudyInFather( Face_6, Edge_16, 'Edge_16' )
geompy.addToStudy( Thickness_1, 'Thickness_1' )
geompy.addToStudy( NoExtraEdges_1, 'NoExtraEdges_1' )
geompy.addToStudyInFather( NoExtraEdges_1, Face_7, 'Face_7' )
geompy.addToStudyInFather( NoExtraEdges_1, Face_8, 'Face_8' )
geompy.addToStudyInFather( NoExtraEdges_1, Face_9, 'Face_9' )
geompy.addToStudyInFather( NoExtraEdges_1, Face_10, 'Face_10' )
geompy.addToStudyInFather( NoExtraEdges_1, Face_11, 'Face_11' )
geompy.addToStudyInFather( NoExtraEdges_1, Face_12, 'Face_12' )
geompy.addToStudyInFather( Face_9, Edge_17, 'Edge_17' )
geompy.addToStudyInFather( Face_9, Edge_18, 'Edge_18' )
geompy.addToStudyInFather( Face_9, Edge_19, 'Edge_19' )
geompy.addToStudyInFather( Face_9, Edge_20, 'Edge_20' )
geompy.addToStudy( Shell_1, 'Shell_1' )
geompy.addToStudyInFather( Shell_1, Edge_21, 'Edge_21' )
geompy.addToStudyInFather( Shell_1, Edge_22, 'Edge_22' )
geompy.addToStudyInFather( Shell_1, Edge_23, 'Edge_23' )
geompy.addToStudyInFather( Shell_1, Edge_24, 'Edge_24' )
geompy.addToStudyInFather( Shell_1, Edge_25, 'Edge_25' )
geompy.addToStudyInFather( Shell_1, Edge_26, 'Edge_26' )
geompy.addToStudyInFather( Shell_1, Edge_27, 'Edge_27' )
geompy.addToStudyInFather( Shell_1, Edge_28, 'Edge_28' )
geompy.addToStudyInFather( Shell_1, Edge_29, 'Edge_29' )
geompy.addToStudyInFather( Shell_1, Edge_30, 'Edge_30' )
geompy.addToStudyInFather( Shell_1, Edge_31, 'Edge_31' )
geompy.addToStudyInFather( Shell_1, Edge_32, 'Edge_32' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Shell_1)
Regular_1D = Mesh_1.Segment()
placeholder = Regular_1D.NumberOfSegments(15)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_1 = Mesh_1.Segment(geom=Edge_21)
nsubs_radial = Regular_1D_1.NumberOfSegments(3)
Propagation_of_1D_Hyp = Regular_1D_1.Propagation()
Regular_1D_2 = Mesh_1.Segment(geom=smeshObj_1)
#status = Mesh_1.AddHypothesis(smeshObj_2,smeshObj_1) ### smeshObj_2 has not been yet created
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,smeshObj_1)
Regular_1D_3 = Mesh_1.Segment(geom=Edge_23)
nsubs_tang_ypx = Regular_1D_3.NumberOfSegments(5)
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_23)
Regular_1D_4 = Mesh_1.Segment(geom=Edge_22)
nsubs_tang_ymx = Regular_1D_4.NumberOfSegments(4)
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_22)
#Mesh_1.GetMesh().RemoveSubMesh( smeshObj_3 ) ### smeshObj_3 has not been yet created
isDone = Mesh_1.Compute()
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.EDGE,SMESH.FT_FreeBorders,SMESH.FT_Undefined,0,SMESH.FT_LogicalNOT)
aCriteria.append(aCriterion)
isDone = Mesh_1.RemoveElements( [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 36, 37, 38, 39, 40, 41, 42, 43, 44 ] )
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_1_0' )
nbAdd = Group_1_0.AddFrom( Mesh_1.GetMesh() )
Sub_mesh_1 = Regular_1D_1.GetSubMesh()
Sub_mesh_2 = Regular_1D_3.GetSubMesh()
Sub_mesh_3 = Regular_1D_4.GetSubMesh()

## some objects were removed
aStudyBuilder = salome.myStudy.NewBuilder()
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_1))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)

## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(nsubs_radial, 'nsubs_radial')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(placeholder, 'placeholder')
smesh.SetName(nsubs_tang_ypx, 'nsubs_tang_ypx')
smesh.SetName(nsubs_tang_ymx, 'nsubs_tang_ymx')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
