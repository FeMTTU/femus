#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.2.1 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/student/software/femus/applications/OptimalControl/boundary_control_inequality/dirichlet/dirichlet_lifting_external/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("y_axis", 2)
notebook.set("x_axis", 2)
notebook.set("x_axis_ext", 1)
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

sk = geompy.Sketcher2D()
sk.addPoint(0.0000000, 0.0000000)
sk.addSegmentAbsolute(10.0000000, 0.0000000)
sk.addSegmentAbsolute(20.0000000, 0.0000000)
sk.addSegmentAbsolute(20.0000000, 20.0000000)
sk.addSegmentAbsolute(0.0000000, 20.0000000)
sk.close()
Face_1 = sk.face([0, 0, 0, 0, 0, 1, 1, 0, 0])
[geomObj_1] = geompy.SubShapes(Face_1, [6])
twoEdges = geompy.CreateGroup(Face_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(twoEdges, [3, 6])
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
O_1 = geompy.MakeVertex(0, 0, 0)
OX_1 = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY_1 = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ_1 = geompy.MakeVectorDXDYDZ(0, 0, 1)
Vertex_1 = geompy.MakeVertex(0, 0, 0)
Vertex_2 = geompy.MakeVertex(1, 0, 0)
Vertex_3 = geompy.MakeVertex(0, 1, 0)
Vertex_4 = geompy.MakeVertex(1, 1, 0)
Vertex_5 = geompy.MakeVertex(1, 0, 0)
Vertex_6 = geompy.MakeVertex(1, 1, 0)
Vertex_7 = geompy.MakeVertex(1.5, 0, 0)
Vertex_8 = geompy.MakeVertex(1.5, 1, 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_4)
Line_3 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_3)
Line_4 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_1)
Line_5 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_7)
Line_6 = geompy.MakeLineTwoPnt(Vertex_8, Vertex_7)
Line_7 = geompy.MakeLineTwoPnt(Vertex_8, Vertex_6)
Line_8 = geompy.MakeLineTwoPnt(Vertex_6, Vertex_5)
Face_1_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_3, Line_4], 1)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_1_1, geompy.ShapeType["EDGE"], True)
Face_2 = geompy.MakeFaceWires([Line_5, Line_6, Line_7, Line_8], 1)
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Face_2, geompy.ShapeType["EDGE"], True)
Vertex_9 = geompy.MakeVertex(1.5, 0.5, 0)
Arc_1 = geompy.MakeArc(Vertex_2, Vertex_9, Vertex_4)
Face_3 = geompy.MakeFaceWires([Line_2, Arc_1], 0)
[Edge_9,Edge_10] = geompy.ExtractShapes(Face_3, geompy.ShapeType["EDGE"], True)
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudyInFather( Face_1, twoEdges, 'twoEdges' )
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( O_1, 'O' )
geompy.addToStudy( OX_1, 'OX' )
geompy.addToStudy( OY_1, 'OY' )
geompy.addToStudy( OZ_1, 'OZ' )
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
geompy.addToStudy( Face_1_1, 'Face_1' )
geompy.addToStudyInFather( Face_1_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_1_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_1_1, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Face_1_1, Edge_4, 'Edge_4' )
geompy.addToStudy( Face_2, 'Face_2' )
geompy.addToStudyInFather( Face_2, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_2, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_2, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Face_2, Edge_8, 'Edge_8' )
geompy.addToStudy( Vertex_9, 'Vertex_9' )
geompy.addToStudy( Arc_1, 'Arc_1' )
geompy.addToStudy( Face_3, 'Face_3' )
geompy.addToStudyInFather( Face_3, Edge_9, 'Edge_9' )
geompy.addToStudyInFather( Face_3, Edge_10, 'Edge_10' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
Mesh_1 = smesh.Mesh(Face_1_1)
Mesh_2 = smesh.Mesh(Face_2)
theNbElems = Mesh_1.Evaluate(Face_1_1)
Regular_1D_2 = Mesh_1.Segment(geom=Edge_1)
Number_of_Segments_1 = Regular_1D_2.NumberOfSegments("y_axis",None,[])
Regular_1D_3 = Mesh_1.Segment(geom=Edge_2)
Number_of_Segments_2 = Regular_1D_3.NumberOfSegments("x_axis",None,[])
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_4 = Mesh_1.Segment()
Regular_1D_5 = Mesh_2.Segment(geom=Edge_5)
Regular_1D_6 = Mesh_2.Segment(geom=Edge_6)
status = Mesh_2.AddHypothesis(Number_of_Segments_1,Edge_5)
Number_of_Segments_3 = Regular_1D_6.NumberOfSegments("x_axis_ext",None,[])
Propagation_of_1D_Hyp = Regular_1D_6.Propagation()
Regular_1D_7 = Mesh_2.Segment()
Quadrangle_2D_1 = Mesh_2.Quadrangle(algo=smeshBuilder.QUADRANGLE)
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_1)
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_2)
isDone = Mesh_1.Compute()
status = Mesh_2.AddHypothesis(Propagation_of_1D_Hyp,Edge_5)
isDone = Mesh_2.Compute()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Mesh_2.ConvertToQuadratic(0, Mesh_2,True)
Mesh_3 = smesh.Concatenate([ Mesh_1.GetMesh(), Mesh_2.GetMesh() ], 1, 1, 1e-05)
isDone = Mesh_1.RemoveElements( [ 5, 6 ] )
Mesh_3.Clear()
isDone = Mesh_3.RemoveElements( [ 3, 4 ] )
Group_13_0 = Mesh_3.CreateEmptyGroup( SMESH.FACE, 'Group_13_0' )
nbAdd = Group_13_0.Add( [ 15, 16 ] )
Group_12_0 = Mesh_3.CreateEmptyGroup( SMESH.FACE, 'Group_12_0' )
nbAdd = Group_12_0.Add( [ 7, 8, 9, 10 ] )
Group_1_0 = Mesh_3.CreateEmptyGroup( SMESH.EDGE, 'Group_1_0' )
nbAdd = Group_1_0.AddFrom( Mesh_3.GetMesh() )
Number_of_Segments_4 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_4.SetNumberOfSegments( 2 )
Number_of_Segments_4.SetReversedEdges( [] )
Number_of_Segments_4.SetObjectEntry( "Face_3" )
NETGEN_1D_2D = smesh.CreateHypothesis('NETGEN_2D', 'libNETGENEngine.so')
Projection_2D = smesh.CreateHypothesis('Projection_2D')
Number_of_Segments_5 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_5.SetNumberOfSegments( 4 )
Number_of_Segments_6 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_6.SetNumberOfSegments( 6 )
Regular_1D_8 = smesh.CreateHypothesis( "Regular_1D" )
Number_of_Segments_7 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_7.SetNumberOfSegments( 4 )
Number_of_Segments_7.SetReversedEdges( [] )
Number_of_Segments_7.SetObjectEntry( "Face_3" )
Mesh_4 = smesh.Mesh(Face_3)
Regular_1D_9 = Mesh_4.Segment()
RadialQuadrangle_1D2D = Mesh_4.Quadrangle(algo=smeshBuilder.RADIAL_QUAD)
Number_of_Segments_8 = Regular_1D_9.NumberOfSegments(2)
status = Mesh_4.RemoveHypothesis(Number_of_Segments_8)
Number_of_Segments_9 = Regular_1D_9.NumberOfSegments(15)
status = Mesh_4.RemoveHypothesis(Number_of_Segments_9)
Number_of_Segments_10 = Regular_1D_9.NumberOfSegments(4)
status = Mesh_4.RemoveHypothesis(RadialQuadrangle_1D2D)
status = Mesh_4.AddHypothesis(Projection_2D)
status = Mesh_4.RemoveHypothesis(Number_of_Segments_10)
Number_of_Segments_11 = Regular_1D_9.NumberOfSegments(10)
isDone = Mesh_4.Compute()
Sub_mesh_1 = Regular_1D_2.GetSubMesh()
Sub_mesh_4 = Regular_1D_3.GetSubMesh()
Sub_mesh_2 = Mesh_2.GetSubMesh( Edge_5, 'Regular_1D' )
Sub_mesh_3 = Mesh_2.GetSubMesh( Edge_6, 'Regular_1D' )


## Set names of Mesh objects
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Projection_2D, 'Projection_2D')
smesh.SetName(NETGEN_1D_2D, 'NETGEN 1D-2D')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(RadialQuadrangle_1D2D.GetAlgorithm(), 'RadialQuadrangle_1D2D')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Number_of_Segments_5, 'Number of Segments_5')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Number_of_Segments_6, 'Number of Segments_6')
smesh.SetName(Mesh_4.GetMesh(), 'Mesh_4')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Mesh_3.GetMesh(), 'Mesh_3')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(Number_of_Segments_10, 'Number of Segments_10')
smesh.SetName(Number_of_Segments_9, 'Number of Segments_9')
smesh.SetName(Number_of_Segments_8, 'Number of Segments_8')
smesh.SetName(Number_of_Segments_7, 'Number of Segments_7')
smesh.SetName(Number_of_Segments_11, 'Number of Segments_11')
smesh.SetName(Group_13_0, 'Group_13_0')
smesh.SetName(Group_12_0, 'Group_12_0')
smesh.SetName(Group_1_0, 'Group_1_0')
