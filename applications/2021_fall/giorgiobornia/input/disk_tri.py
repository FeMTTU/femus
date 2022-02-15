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

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Disk_1)
Regular_1D = Mesh_1.Segment()
Max_Size_1 = Regular_1D.MaxSize(0.282843)
MEFISTO_2D = Mesh_1.Triangle(algo=smeshBuilder.MEFISTO)
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0)
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_1' )
nbAdd = Group_1_0.AddFrom( Mesh_1.GetMesh() )
Group_1_0.SetName( 'Group_1_0' )
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/gbornia/software/femus/applications/OptimalControl/00_laplacian_all_dims/input/disk_tri.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Mesh_1_rotated = Mesh_1.RotateObjectMakeMesh( Mesh_1, SMESH.AxisStruct( 0, 0, 0, 1, 0, 0 ), 0.785398, 1, 'Mesh_1_rotated' )
[ Group_1_0_1 ] = Mesh_1_rotated.GetGroups()
Mesh_1_rotated_rotated = Mesh_1_rotated.RotateObjectMakeMesh( Mesh_1_rotated, SMESH.AxisStruct( 0, 0, 0, 1, 0, 0 ), 0.785398, 1, 'Mesh_1_rotated_rotated' )
[ Group_1_0_2 ] = Mesh_1_rotated_rotated.GetGroups()
[ Group_1_0_1 ] = Mesh_1_rotated.GetGroups()
[ Group_1_0_2 ] = Mesh_1_rotated_rotated.GetGroups()
smesh.SetName(Mesh_1_rotated, 'Mesh_1_rotated')
try:
  Mesh_1_rotated.ExportMED(r'/home/gbornia/software/femus/applications/OptimalControl/00_laplacian_all_dims/input/disk_tri_45x.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
smesh.SetName(Mesh_1_rotated_rotated, 'Mesh_1_rotated_rotated')
try:
  Mesh_1_rotated_rotated.ExportMED(r'/home/gbornia/software/femus/applications/OptimalControl/00_laplacian_all_dims/input/disk_tri_90x.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(Group_1_0_1, 'Group_1_0')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(MEFISTO_2D.GetAlgorithm(), 'MEFISTO_2D')
smesh.SetName(Max_Size_1, 'Max Size_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Mesh_1_rotated_rotated.GetMesh(), 'Mesh_1_rotated_rotated')
smesh.SetName(Mesh_1_rotated.GetMesh(), 'Mesh_1_rotated')
smesh.SetName(Group_1_0_2, 'Group_1_0')
smesh.SetName(Group_1_0, 'Group_1_0')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
