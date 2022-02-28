#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.3.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/dn27/users_Linux/eap/salome/tmp')

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
Divided_Disk_1 = geompy.MakeDividedDisk(1, 1, GEOM.SQUARE)
[Edge_1,Edge_2,Edge_3,Edge_4,Edge_5,Edge_6,Edge_7,Edge_8,Edge_9,Edge_10,Edge_11,Edge_12] = geompy.ExtractShapes(Divided_Disk_1, geompy.ShapeType["EDGE"], True)
Auto_group_for_Sub_mesh_1 = geompy.CreateGroup(Divided_Disk_1, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Sub_mesh_1, [Edge_2, Edge_3])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Divided_Disk_1, 'Divided Disk_1' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_9, 'Edge_9' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_10, 'Edge_10' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_11, 'Edge_11' )
geompy.addToStudyInFather( Divided_Disk_1, Edge_12, 'Edge_12' )
geompy.addToStudyInFather( Divided_Disk_1, Auto_group_for_Sub_mesh_1, 'Auto_group_for_Sub-mesh_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Divided_Disk_1)
Regular_1D = Mesh_1.Segment()
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_1 = Mesh_1.Segment(geom=Auto_group_for_Sub_mesh_1)
central_dir_ymx = Regular_1D_1.NumberOfSegments(4)
Regular_1D_2 = Mesh_1.Segment(geom=Edge_1)
nsubs_radial = Regular_1D_2.NumberOfSegments(3)
Propagation_of_1D_Hyp = Regular_1D_2.Propagation()
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Auto_group_for_Sub_mesh_1)
isDone = Mesh_1.Compute()
Sub_mesh_1 = Regular_1D_1.GetSubMesh()
Sub_mesh_2 = Regular_1D_2.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(nsubs_radial, 'nsubs_radial')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(central_dir_ymx, 'central_dir_ymx')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
