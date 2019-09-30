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
Divided_Disk_1 = geompy.MakeDividedDisk(1, 1, GEOM.SQUARE)
[Face_1,Face_2,Face_3,Face_4,Face_5] = geompy.ExtractShapes(Divided_Disk_1, geompy.ShapeType["FACE"], True)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Face_2, geompy.ShapeType["EDGE"], True)
[Edge_9,Edge_10,Edge_11,Edge_12] = geompy.ExtractShapes(Face_3, geompy.ShapeType["EDGE"], True)
[Edge_13,Edge_14,Edge_15,Edge_16] = geompy.ExtractShapes(Face_4, geompy.ShapeType["EDGE"], True)
[Edge_17,Edge_18,Edge_19,Edge_20] = geompy.ExtractShapes(Face_5, geompy.ShapeType["EDGE"], True)
geomObj_1 = geompy.MakeGlueEdges([Face_1, Face_2], 1e-07)
geomObj_2 = geompy.MakeGlueEdges([Face_1, Face_2], 1e-07)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Divided_Disk_1, 'Divided Disk_1' )
geompy.addToStudyInFather( Divided_Disk_1, Face_1, 'Face_1' )
geompy.addToStudyInFather( Divided_Disk_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Divided_Disk_1, Face_3, 'Face_3' )
geompy.addToStudyInFather( Divided_Disk_1, Face_4, 'Face_4' )
geompy.addToStudyInFather( Divided_Disk_1, Face_5, 'Face_5' )
geompy.addToStudyInFather( Face_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_1, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Face_1, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( Face_2, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_2, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_2, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Face_2, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Face_3, Edge_9, 'Edge_9' )
geompy.addToStudyInFather( Face_3, Edge_10, 'Edge_10' )
geompy.addToStudyInFather( Face_3, Edge_11, 'Edge_11' )
geompy.addToStudyInFather( Face_3, Edge_12, 'Edge_12' )
geompy.addToStudyInFather( Face_4, Edge_13, 'Edge_13' )
geompy.addToStudyInFather( Face_4, Edge_14, 'Edge_14' )
geompy.addToStudyInFather( Face_4, Edge_15, 'Edge_15' )
geompy.addToStudyInFather( Face_4, Edge_16, 'Edge_16' )
geompy.addToStudyInFather( Face_5, Edge_17, 'Edge_17' )
geompy.addToStudyInFather( Face_5, Edge_18, 'Edge_18' )
geompy.addToStudyInFather( Face_5, Edge_19, 'Edge_19' )
geompy.addToStudyInFather( Face_5, Edge_20, 'Edge_20' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Face_1)
Regular_1D = Mesh_1.Segment()
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_1 = Mesh_1.Segment(geom=Edge_1)
nsubs_radial = Regular_1D_1.NumberOfSegments(3)
Regular_1D_2 = Mesh_1.Segment(geom=Edge_3)
central_dir_ymx = Regular_1D_2.NumberOfSegments(4)
placeholder = Regular_1D.NumberOfSegments(15)
Propagation_of_1D_Hyp = smesh.CreateHypothesis('Propagation')
Propagation_of_1D_Hyp_1 = Regular_1D_2.Propagation()
isDone = Mesh_1.Compute()
Sub_mesh_1 = Regular_1D_1.GetSubMesh()
Sub_mesh_2 = Regular_1D_2.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(central_dir_ymx, 'central_dir_ymx')
smesh.SetName(placeholder, 'placeholder')
smesh.SetName(nsubs_radial, 'nsubs_radial')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Propagation_of_1D_Hyp_1, 'Propagation of 1D Hyp. on Opposite Edges_2')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
