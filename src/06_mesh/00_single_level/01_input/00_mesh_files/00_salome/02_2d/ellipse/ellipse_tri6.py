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


geompy = geomBuilder.New()

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

smesh = smeshBuilder.New()

Max_Size_2 = smesh.CreateHypothesis('MaxLength')
Max_Size_2.SetLength( 0.447214 )
Mesh_3 = smesh.Mesh(Face_1)
status = Mesh_3.AddHypothesis(Max_Size_2)
status = Mesh_3.AddHypothesis(Regular_1D)
status = Mesh_3.AddHypothesis(NETGEN_2D)
isDone = Mesh_3.Compute()
Mesh_3.ConvertToQuadratic(0)


Group_1_0 = Mesh_3.CreateEmptyGroup( SMESH.EDGE, 'Group_1' )
nbAdd = Group_1_0.AddFrom( Mesh_3.GetMesh() )
Group_1_0.SetName( 'Group_1_0' )


## Set names of Mesh objects
smesh.SetName(Max_Size_2, 'Max Size_2')
smesh.SetName(Mesh_3.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')

