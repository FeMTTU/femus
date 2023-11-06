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
notebook.set("L", 5)
notebook.set("radius", 1)
notebook.set("n_r", 7)
notebook.set("n_s", 10)
notebook.set("L", 5)
notebook.set("radius", 1)
notebook.set("n_r", 7)
notebook.set("n_s", 10)
notebook.set("L", 5)
notebook.set("radius", 1)
notebook.set("n_r", 7)
notebook.set("n_s", 10)
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
Vertex_1 = geompy.MakeVertex(0, 0, 0)
Vertex_2 = geompy.MakeVertex("radius", 0, 0)
Vertex_3 = geompy.MakeVertex(0, "radius", 0)
Arc_1 = geompy.MakeArcCenter(Vertex_1, Vertex_3, Vertex_2,False)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_3)
Line_1_vertex_2 = geompy.GetSubShape(Line_1, [2])
Line_2 = geompy.MakeLineTwoPnt(Line_1_vertex_2, Vertex_2)
Face_1 = geompy.MakeFaceWires([Arc_1, Line_1, Line_2], 1)
[Edge_1,Edge_2,Edge_3] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)
[Edge_1, Edge_2, Edge_3] = geompy.GetExistingSubObjects(Face_1, False)
[Edge_1, Edge_2, Edge_3] = geompy.GetExistingSubObjects(Face_1, False)
[Edge_1, Edge_2, Edge_3] = geompy.GetExistingSubObjects(Face_1, False)
[Edge_1, Edge_2, Edge_3] = geompy.GetExistingSubObjects(Face_1, False)
[Edge_1, Edge_2, Edge_3] = geompy.GetExistingSubObjects(Face_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Arc_1, 'Arc_1' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudyInFather( Line_1, Line_1_vertex_2, 'Line_1:vertex_2' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudyInFather( Face_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_1, Edge_3, 'Edge_3' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

NETGEN_2D_Parameters_1 = smesh.CreateHypothesisByAverageLength( 'NETGEN_Parameters_2D', 'NETGENEngine', 0.282843, 0 )
NETGEN_1D_2D = smesh.CreateHypothesis('NETGEN_2D', 'libNETGENEngine.so')
Regular_1D = smesh.CreateHypothesis('Regular_1D')
Number_of_Segments_1 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_1.SetNumberOfSegments( "n_r" )
Number_of_Segments_1.SetReversedEdges( [] )
Number_of_Segments_1.SetObjectEntry( "Face_1" )
Number_of_Segments_2 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_2.SetNumberOfSegments( "n_s" )
Number_of_Segments_2.SetReversedEdges( [] )
Number_of_Segments_2.SetObjectEntry( "Face_1" )
#hyp_14.SetLength( 0.141421 ) ### not created Object
NETGEN_2D_Parameters_2 = smesh.CreateHypothesisByAverageLength( 'NETGEN_Parameters_2D', 'NETGENEngine', 0.141421, 0 )
try:
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Number_of_Segments_1.SetNumberOfSegments( "n_r" )
Propagation_of_1D_Hyp = smesh.CreateHypothesis('Propagation')
Mesh_1 = smesh.Mesh(Face_1)
Quadratic_Mesh_1 = smesh.CreateHypothesis('QuadraticMesh')
status = Mesh_1.AddHypothesis(Regular_1D,Edge_1)
status = Mesh_1.AddHypothesis(Number_of_Segments_1,Edge_1)
status = Mesh_1.AddHypothesis(Regular_1D,Edge_2)
status = Mesh_1.AddHypothesis(Number_of_Segments_1,Edge_2)
status = Mesh_1.AddHypothesis(Regular_1D,Edge_3)
status = Mesh_1.AddHypothesis(Number_of_Segments_1,Edge_3)
#hyp_17.SetLength( 0.141421 ) ### not created Object
NETGEN_2D_Parameters_3 = smesh.CreateHypothesisByAverageLength( 'NETGEN_Parameters_2D', 'NETGENEngine', 0.141421, 0 )
status = Mesh_1.AddHypothesis(NETGEN_1D_2D)
status = Mesh_1.AddHypothesis( Face_1, NETGEN_2D_Parameters_3 )
status = Mesh_1.RemoveHypothesis(NETGEN_2D_Parameters_3)
status = Mesh_1.AddHypothesis( Face_1, NETGEN_2D_Parameters_1 )
Group_1_0 = Mesh_1.GroupOnGeom(Edge_1,'Group_1_0',SMESH.EDGE)
[ smeshObj_1, smeshObj_2, smeshObj_3, Group_1_0 ] = Mesh_1.GetGroups()
Group_2_0 = Mesh_1.GroupOnGeom(Edge_2,'Group_2_0',SMESH.EDGE)
[ smeshObj_1, smeshObj_2, smeshObj_3, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Group_3_0 = Mesh_1.GroupOnGeom(Edge_3,'Group_3_0',SMESH.EDGE)
[ smeshObj_1, smeshObj_2, smeshObj_3, Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/max/software/femus/applications/2021_fall/shaflynn/input/Mesh_1_assignment_1_triangular.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/max/software/femus/applications/2021_fall/shaflynn/input/Mesh_1_assignment_1_triangular.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
isDone = Mesh_1.Compute()
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0)
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
Sub_mesh_1 = Mesh_1.GetSubMesh( Edge_1, 'Sub-mesh_1' )
Sub_mesh_2 = Mesh_1.GetSubMesh( Edge_2, 'Sub-mesh_2' )
Sub_mesh_3 = Mesh_1.GetSubMesh( Edge_3, 'Sub-mesh_3' )

## some objects were removed
aStudyBuilder = salome.myStudy.NewBuilder()
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_2))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_3))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_1))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)

## Set names of Mesh objects
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(NETGEN_1D_2D, 'NETGEN 1D-2D')
smesh.SetName(Regular_1D, 'Regular_1D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(NETGEN_2D_Parameters_3, 'NETGEN 2D Parameters_3')
smesh.SetName(Quadratic_Mesh_1, 'Quadratic Mesh_1')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(NETGEN_2D_Parameters_2, 'NETGEN 2D Parameters_2')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
