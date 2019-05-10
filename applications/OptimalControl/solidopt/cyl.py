#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.2.1 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/applications/OptimalControl/solidopt')

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
Cylinder_1 = geompy.MakeCylinder(O, OZ, 1, 1)
[Face_1,Face_2,Face_3] = geompy.ExtractShapes(Cylinder_1, geompy.ShapeType["FACE"], True)
[Edge_1,Edge_2,Edge_3] = geompy.ExtractShapes(Cylinder_1, geompy.ShapeType["EDGE"], True)
Face_4 = geompy.MakeFaceHW(0.8, 0.8, 1)
[Vertex_1,Vertex_2,Vertex_3,Vertex_4] = geompy.ExtractShapes(Face_4, geompy.ShapeType["VERTEX"], True)
[Edge_4,Edge_5,Edge_6,Edge_7] = geompy.ExtractShapes(Face_4, geompy.ShapeType["EDGE"], True)
Line_1 = geompy.MakeLineTwoPnt(O, Vertex_1)
ExtendedEdge_1 = geompy.ExtendEdge(Line_1, 0, 3)
Line_2 = geompy.MakeLineTwoPnt(O, Vertex_2)
Line_3 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_4 = geompy.MakeLineTwoPnt(O, Vertex_4)
ExtendedEdge_2 = geompy.ExtendEdge(Line_2, 0, 2)
ExtendedEdge_3 = geompy.ExtendEdge(Line_3, 0, 2)
ExtendedEdge_4 = geompy.ExtendEdge(Line_4, 0, 2)
Intersection_1 = geompy.MakeSection(Edge_1, ExtendedEdge_4, True)
Intersection_2 = geompy.MakeSection(Edge_1, ExtendedEdge_1, True)
Intersection_3 = geompy.MakeSection(ExtendedEdge_3, Edge_1, True)
Intersection_4 = geompy.MakeSection(ExtendedEdge_2, Edge_1, True)
Arc_1 = geompy.MakeArcCenter(O, Intersection_1, Intersection_4,False)
Arc_2 = geompy.MakeArcCenter(O, Intersection_4, Intersection_2,False)
Arc_3 = geompy.MakeArcCenter(O, Intersection_2, Intersection_3,False)
Arc_4 = geompy.MakeArcCenter(O, Intersection_3, Intersection_1,False)
Arc_2_vertex_3 = geompy.GetSubShape(Arc_2, [3])
Line_5 = geompy.MakeLineTwoPnt(Vertex_1, Arc_2_vertex_3)
Arc_2_vertex_2 = geompy.GetSubShape(Arc_2, [2])
Line_6 = geompy.MakeLineTwoPnt(Arc_2_vertex_2, Vertex_2)
Arc_4_vertex_2 = geompy.GetSubShape(Arc_4, [2])
Line_7 = geompy.MakeLineTwoPnt(Arc_4_vertex_2, Vertex_3)
Arc_1_vertex_2 = geompy.GetSubShape(Arc_1, [2])
Line_8 = geompy.MakeLineTwoPnt(Arc_1_vertex_2, Vertex_4)
Face_5 = geompy.MakeFaceWires([Arc_4, Line_7, Line_8, Edge_7], 1)
[Edge_8,Edge_9,Edge_10,Edge_11] = geompy.ExtractShapes(Face_5, geompy.ShapeType["EDGE"], True)
Face_6 = geompy.MakeFaceWires([Arc_1, Line_6, Edge_6, Edge_10], 1)
[Edge_12,Edge_13,Edge_14,Edge_15] = geompy.ExtractShapes(Face_6, geompy.ShapeType["EDGE"], True)
Face_7 = geompy.MakeFaceWires([Arc_3, Line_5, Edge_5, Edge_9], 1)
[Edge_16,Edge_17,Edge_18,Edge_19] = geompy.ExtractShapes(Face_7, geompy.ShapeType["EDGE"], True)
Face_8 = geompy.MakeFaceWires([Arc_2, Edge_4, Edge_12, Edge_16], 1)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudyInFather( Cylinder_1, Face_1, 'Face_1' )
geompy.addToStudyInFather( Cylinder_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Cylinder_1, Face_3, 'Face_3' )
geompy.addToStudy( Face_4, 'Face_4' )
geompy.addToStudyInFather( Face_4, Vertex_1, 'Vertex_1' )
geompy.addToStudyInFather( Face_4, Vertex_2, 'Vertex_2' )
geompy.addToStudyInFather( Face_4, Vertex_3, 'Vertex_3' )
geompy.addToStudyInFather( Face_4, Vertex_4, 'Vertex_4' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( ExtendedEdge_1, 'ExtendedEdge_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudyInFather( Cylinder_1, Edge_1, 'Edge_1' )
geompy.addToStudy( ExtendedEdge_2, 'ExtendedEdge_2' )
geompy.addToStudy( ExtendedEdge_3, 'ExtendedEdge_3' )
geompy.addToStudy( ExtendedEdge_4, 'ExtendedEdge_4' )
geompy.addToStudyInFather( Cylinder_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Cylinder_1, Edge_3, 'Edge_3' )
geompy.addToStudy( Intersection_1, 'Intersection_1' )
geompy.addToStudy( Intersection_2, 'Intersection_2' )
geompy.addToStudy( Intersection_3, 'Intersection_3' )
geompy.addToStudy( Intersection_4, 'Intersection_4' )
geompy.addToStudy( Arc_2, 'Arc_2' )
geompy.addToStudy( Arc_1, 'Arc_1' )
geompy.addToStudy( Arc_3, 'Arc_3' )
geompy.addToStudy( Arc_4, 'Arc_4' )
geompy.addToStudyInFather( Arc_2, Arc_2_vertex_3, 'Arc_2:vertex_3' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudyInFather( Arc_2, Arc_2_vertex_2, 'Arc_2:vertex_2' )
geompy.addToStudy( Line_6, 'Line_6' )
geompy.addToStudyInFather( Arc_4, Arc_4_vertex_2, 'Arc_4:vertex_2' )
geompy.addToStudy( Line_7, 'Line_7' )
geompy.addToStudyInFather( Arc_1, Arc_1_vertex_2, 'Arc_1:vertex_2' )
geompy.addToStudy( Line_8, 'Line_8' )
geompy.addToStudyInFather( Face_4, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( Face_4, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_4, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_4, Edge_7, 'Edge_7' )
geompy.addToStudy( Face_5, 'Face_5' )
geompy.addToStudyInFather( Face_5, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Face_5, Edge_9, 'Edge_9' )
geompy.addToStudyInFather( Face_5, Edge_10, 'Edge_10' )
geompy.addToStudyInFather( Face_5, Edge_11, 'Edge_11' )
geompy.addToStudy( Face_6, 'Face_6' )
geompy.addToStudy( Face_7, 'Face_7' )
geompy.addToStudyInFather( Face_6, Edge_12, 'Edge_12' )
geompy.addToStudyInFather( Face_6, Edge_13, 'Edge_13' )
geompy.addToStudyInFather( Face_6, Edge_14, 'Edge_14' )
geompy.addToStudyInFather( Face_6, Edge_15, 'Edge_15' )
geompy.addToStudyInFather( Face_7, Edge_16, 'Edge_16' )
geompy.addToStudyInFather( Face_7, Edge_17, 'Edge_17' )
geompy.addToStudyInFather( Face_7, Edge_18, 'Edge_18' )
geompy.addToStudyInFather( Face_7, Edge_19, 'Edge_19' )
geompy.addToStudy( Face_8, 'Face_8' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
Mesh_1 = smesh.Mesh(Cylinder_1)
Number_of_Segments_1 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_1.SetNumberOfSegments( 2 )
NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_2D,geom=Face_1)
status = Mesh_1.RemoveHypothesis(NETGEN_2D,Face_1)
Regular_1D = Mesh_1.Segment(geom=Face_1)
Max_Size_1 = Regular_1D.MaxSize(0.3)
MEFISTO_2D = Mesh_1.Triangle(algo=smeshBuilder.MEFISTO,geom=Face_1)
Prism_3D = Mesh_1.Prism()
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_1 = Mesh_1.Segment()
Number_of_Segments_2 = Regular_1D_1.NumberOfSegments(3)
Number_of_Segments_2.SetNumberOfSegments( 1 )
isDone = Mesh_1.Compute()
Sub_mesh_1 = NETGEN_2D.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(MEFISTO_2D.GetAlgorithm(), 'MEFISTO_2D')
smesh.SetName(Max_Size_1, 'Max Size_1')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Prism_3D.GetAlgorithm(), 'Prism_3D')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
