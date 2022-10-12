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
sys.path.insert( 0, r'/home/gbornia/software/femus')

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
Disk_1 = geompy.MakeDiskR(1, 1)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Disk_1, 'Disk_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
Number_of_Layers_1 = smesh.CreateHypothesis('NumberOfLayers2D')
Number_of_Layers_1.SetNumberOfLayers( 4 )
RadialQuadrangle_1D2D = smesh.CreateHypothesis('RadialQuadrangle_1D2D')
Number_of_Layers_2 = smesh.CreateHypothesis('NumberOfLayers2D')
Number_of_Layers_2.SetNumberOfLayers( 2 )
Number_of_Layers_3 = smesh.CreateHypothesis('NumberOfLayers2D')
Number_of_Layers_3.SetNumberOfLayers( 2 )
Mesh_1 = smesh.Mesh(Disk_1)
NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_2D)
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0)
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_1' )
nbAdd = Group_1_0.AddFrom( Mesh_1.GetMesh() )
Group_1_0.SetName( 'Group_1_0' )
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED( r'/home/gbornia/software/femus/circle_tri6.med', 0, SMESH.MED_V2_2, 1, None ,0)
  pass
except:
  print 'ExportToMEDX() failed. Invalid file name?'


## Set names of Mesh objects
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(RadialQuadrangle_1D2D, 'RadialQuadrangle_1D2D')
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(Number_of_Layers_2, 'Number of Layers_2')
smesh.SetName(Number_of_Layers_3, 'Number of Layers_3')
smesh.SetName(Number_of_Layers_1, 'Number of Layers_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
