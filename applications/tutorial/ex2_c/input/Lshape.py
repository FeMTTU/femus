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
sys.path.insert( 0, r'/home/gbornia/software/femus/applications/tutorial/ex_time/input')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

length_x = 1.
length_y = 5.

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Vertex_1 = geompy.MakeVertex(0, 0, 0)
Vertex_2 = geompy.MakeVertex(1, 0, 0)
Vertex_3 = geompy.MakeVertex(0, length_y, 0)
Vertex_4 = geompy.MakeVertex(1, 0.5, 0)
Vertex_5 = geompy.MakeVertex(0.5, 0.5, 0)
Vertex_6 = geompy.MakeVertex(0.5, length_y, 0)
Vertex_7 = geompy.MakeVertex(0, 0.5, 0)
Line_2 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_6)
Line_3 = geompy.MakeLineTwoPnt(Vertex_6, Vertex_5)
Line_4 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_4)
Line_5 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_2)
Line_6 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_1)
Line_7 = geompy.MakeLineTwoPnt(Vertex_7, Vertex_5)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_7)
Line_8 = geompy.MakeLineTwoPnt(Vertex_7, Vertex_3)
Face_1 = geompy.MakeFaceWires([Line_2, Line_3, Line_4, Line_5, Line_6, Line_1, Line_8], 1)
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
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudy( Line_6, 'Line_6' )
geompy.addToStudy( Line_7, 'Line_7' )
geompy.addToStudy( Line_8, 'Line_8' )
geompy.addToStudy( Face_1, 'Face_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
Number_of_Segments_1 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_1.SetNumberOfSegments( 2 )
Regular_1D = smesh.CreateHypothesis('Regular_1D')
Quadrangle_2D = smesh.CreateHypothesis('Quadrangle_2D')
Mesh_1 = smesh.Mesh(Face_1)
NETGEN_1D_2D = smesh.CreateHypothesis('NETGEN_2D', 'NETGENEngine')
Max_Size_1 = smesh.CreateHypothesis('MaxLength')
Max_Size_1.SetLength( 0.141421 )
status = Mesh_1.AddHypothesis(Regular_1D)
status = Mesh_1.AddHypothesis(Max_Size_1)
MEFISTO_2D = Mesh_1.Triangle(algo=smeshBuilder.MEFISTO)
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0)
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_1_0' )
nbAdd = Group_1_0.AddFrom( Mesh_1.GetMesh() )
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED( r'/home/gbornia/software/femus/applications/tutorial/ex_time/input/Lshape.med', 0, SMESH.MED_V2_2, 1, None ,0)
  pass
except:
  print 'ExportToMEDX() failed. Invalid file name?'


## Set names of Mesh objects
smesh.SetName(Regular_1D, 'Regular_1D')
smesh.SetName(NETGEN_1D_2D, 'NETGEN 1D-2D')
smesh.SetName(Quadrangle_2D, 'Quadrangle_2D')
smesh.SetName(Max_Size_1, 'Max Size_1')
smesh.SetName(MEFISTO_2D.GetAlgorithm(), 'MEFISTO_2D')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
