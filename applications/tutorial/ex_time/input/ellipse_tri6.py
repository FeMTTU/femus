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

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Ellipse_1 = geompy.MakeEllipse(None, None, 2, 1)
Face_1 = geompy.MakeFaceWires([Ellipse_1], 1)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Ellipse_1, 'Ellipse_1' )
geompy.addToStudy( Face_1, 'Face_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh(Face_1)
NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_2D)
Length_From_Edges_1 = smesh.CreateHypothesis('LengthFromEdges')
status = Mesh_1.RemoveHypothesis(NETGEN_2D)
Regular_1D = Mesh_1.Segment()
Max_Size_1 = Regular_1D.MaxSize(0.447214)
MEFISTO_2D = Mesh_1.Triangle(algo=smeshBuilder.MEFISTO)
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0)
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_1' )
nbAdd = Group_1_0.AddFrom( Mesh_1.GetMesh() )
Group_1_0.SetName( 'Group_1_0' )
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED( r'/home/gbornia/software/femus/applications/tutorial/ex_time/input/ellipse_tri6.med', 0, SMESH.MED_V2_2, 1, None ,0)
  pass
except:
  print 'ExportToMEDX() failed. Invalid file name?'


## Set names of Mesh objects
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(MEFISTO_2D.GetAlgorithm(), 'MEFISTO_2D')
smesh.SetName(Max_Size_1, 'Max Size_1')
smesh.SetName(Length_From_Edges_1, 'Length From Edges_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
