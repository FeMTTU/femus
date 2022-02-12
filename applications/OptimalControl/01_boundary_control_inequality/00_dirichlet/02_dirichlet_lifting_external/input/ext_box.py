# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.5.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, r'/home/student/software/femus/applications/OptimalControl/boundary_control_inequality/dirichlet/dirichlet_lifting_external/input')

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
Vertex_2 = geompy.MakeVertex(1, 0, 0)
Vertex_3 = geompy.MakeVertex(0, 1, 0)
Vertex_4 = geompy.MakeVertex(1, 1, 0)
Vertex_5 = geompy.MakeVertex(1, 0.25, 0)
Vertex_6 = geompy.MakeVertex(1, 0.75, 0)
Vertex_7 = geompy.MakeVertex(1.25, 0.25, 0)
Vertex_8 = geompy.MakeVertex(1.25, 0.75, 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_4)
Line_3 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_3)
Line_4 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_1)
Line_5 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_7)
Line_6 = geompy.MakeLineTwoPnt(Vertex_8, Vertex_7)
Line_7 = geompy.MakeLineTwoPnt(Vertex_8, Vertex_6)
Line_8 = geompy.MakeLineTwoPnt(Vertex_6, Vertex_5)
Face_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_3, Line_4], 1)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)
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
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Vertex_7, 'Vertex_7' )
geompy.addToStudy( Vertex_8, 'Vertex_8' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudy( Line_6, 'Line_6' )
geompy.addToStudy( Line_7, 'Line_7' )
geompy.addToStudy( Line_8, 'Line_8' )
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

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
Mesh_1 = smesh.Mesh(Face_1)
Regular_1D = Mesh_1.Segment(geom=Edge_1)
Number_of_Segments_1 = Regular_1D.NumberOfSegments(4,None,[])
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Number_of_Segments_2 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_2.SetReversedEdges( [] )
Number_of_Segments_2.SetObjectEntry( "Face_1" )
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_1 = Mesh_1.Segment()
Number_of_Segments_3 = Regular_1D_1.NumberOfSegments(4,None,[])
Number_of_Segments_2.SetNumberOfSegments( 4 )
isDone = Mesh_1.Compute()
Mesh_2 = smesh.Mesh(Face_2)
Regular_1D_2 = Mesh_2.Segment(geom=Edge_5)
Number_of_Segments_4 = Regular_1D_2.NumberOfSegments(2,None,[])
Regular_1D_3 = Mesh_2.Segment()
Number_of_Segments_5 = Regular_1D_3.NumberOfSegments(1,None,[])
Quadrangle_2D_1 = Mesh_2.Quadrangle(algo=smeshBuilder.QUADRANGLE)
status = Mesh_2.AddHypothesis(Propagation_of_1D_Hyp,Edge_5)
isDone = Mesh_2.Compute()
#try:
  #pass
#except:
  #print 'ExportToMEDX() failed. Invalid file name?'
#try:
  #pass
#except:
  #print 'ExportToMEDX() failed. Invalid file name?'
#try:
  #pass
#except:
  #print 'ExportToMEDX() failed. Invalid file name?'
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Mesh_2.ConvertToQuadratic(0, Mesh_2,True)
isDone = Mesh_1.Compute()
isDone = Mesh_2.Compute()
Mesh_4 = smesh.Concatenate([ Mesh_1.GetMesh(), Mesh_2.GetMesh() ], 1, 1, 1e-05)
isDone = Mesh_4.RemoveElements( [ 6, 7 ] )
Group_12_0 = Mesh_4.CreateEmptyGroup( SMESH.FACE, 'Group_12_0' )
nbAdd = Group_12_0.Add( [ 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 ] )
Group_13_0 = Mesh_4.CreateEmptyGroup( SMESH.FACE, 'Group_13_0' )
nbAdd = Group_13_0.Add( [ 35, 36 ] )
Group_10_0 = Mesh_4.CreateEmptyGroup( SMESH.EDGE, 'Group_1' )
nbAdd = Group_10_0.AddFrom( Mesh_4.GetMesh() )
Group_10_0.SetName( 'Group_10_0' )
Sub_mesh_1 = Regular_1D.GetSubMesh()
Regular_1D_4 = Mesh_1.GetSubMesh( Edge_2, 'Regular_1D' )
Sub_mesh_2 = Regular_1D_2.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_5, 'Number of Segments_5')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')
smesh.SetName(Mesh_4.GetMesh(), 'Mesh_4')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Group_12_0, 'Group_12_0')
smesh.SetName(Group_13_0, 'Group_13_0')
smesh.SetName(Group_10_0, 'Group_10_0')
smesh.SetName(Regular_1D_4, 'Regular_1D')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
