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
Ellipse_2 = geompy.MakeEllipse(None, None, 0.5, 0.25, OY)
Face_2 = geompy.MakeFaceWires([Ellipse_2], 1)
Cut_1 = geompy.MakeCutList(Face_1, [Face_2], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Ellipse_1, 'Ellipse_1' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Ellipse_2, 'Ellipse_2' )
geompy.addToStudy( Face_2, 'Face_2' )
geompy.addToStudy( Cut_1, 'Cut_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
NETGEN_2D = smesh.CreateHypothesis('NETGEN_2D_ONLY', 'NETGENEngine')
Length_From_Edges_1 = smesh.CreateHypothesis('LengthFromEdges')
Max_Size_1 = smesh.CreateHypothesis('MaxLength')
Max_Size_1.SetLength( 0.447214 )
MEFISTO_2D = smesh.CreateHypothesis('MEFISTO_2D')
Regular_1D = smesh.CreateHypothesis('Regular_1D')
try:
  pass
except:
  print 'ExportToMEDX() failed. Invalid file name?'
RadialQuadrangle_1D2D = smesh.CreateHypothesis('RadialQuadrangle_1D2D')
Max_Size_2 = smesh.CreateHypothesis('MaxLength')
Max_Size_2.SetLength( 1 )
Max_Size_3 = smesh.CreateHypothesis('MaxLength')
Max_Size_3.SetLength( 0.447214 )
Max_Size_4 = smesh.CreateHypothesis('MaxLength')
Max_Size_4.SetLength( 0.447214 )
Mesh_1 = smesh.Mesh(Cut_1)
status = Mesh_1.AddHypothesis(Max_Size_4)
status = Mesh_1.AddHypothesis(Regular_1D)
status = Mesh_1.AddHypothesis(MEFISTO_2D)
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0)
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_1' )
nbAdd = Group_1_0.Add( [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 ] )
Group_1_0.SetName( 'Group_1_0' )
Group_2_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_2_0' )
nbAdd = Group_2_0.Add( [ 23, 24, 25, 26, 27, 28 ] )
[ Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED( r'/home/gbornia/software/femus/applications/tutorial/ex_time/input/ellipse_with_hole_tri6.med', 0, SMESH.MED_V2_2, 1, None ,0)
  pass
except:
  print 'ExportToMEDX() failed. Invalid file name?'


## Set names of Mesh objects
smesh.SetName(NETGEN_2D, 'NETGEN 2D')
smesh.SetName(Regular_1D, 'Regular_1D')
smesh.SetName(MEFISTO_2D, 'MEFISTO_2D')
smesh.SetName(Max_Size_1, 'Max Size_1')
smesh.SetName(Max_Size_2, 'Max Size_2')
smesh.SetName(RadialQuadrangle_1D2D, 'RadialQuadrangle_1D2D')
smesh.SetName(Length_From_Edges_1, 'Length From Edges_1')
smesh.SetName(Max_Size_3, 'Max Size_3')
smesh.SetName(Max_Size_4, 'Max Size_4')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
