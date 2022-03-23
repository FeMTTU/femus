#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.8.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/unittests/test_mesh_read_write/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("dx", 1)
notebook.set("dy", 2)
notebook.set("dz", 3)
notebook.set("mesh_nx", 2)
notebook.set("mesh_ny", 3)
notebook.set("mesh_nz", 5)
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
Box_1 = geompy.MakeBoxDXDYDZ("dx", "dy", "dz")
[Face_1,Face_2,Face_3,Face_4,Face_5,Face_6] = geompy.ExtractShapes(Box_1, geompy.ShapeType["FACE"], True)
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Face_2, geompy.ShapeType["EDGE"], True)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_3, geompy.ShapeType["EDGE"], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudyInFather( Box_1, Face_1, 'Face_1' )
geompy.addToStudyInFather( Box_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Box_1, Face_3, 'Face_3' )
geompy.addToStudyInFather( Box_1, Face_4, 'Face_4' )
geompy.addToStudyInFather( Box_1, Face_5, 'Face_5' )
geompy.addToStudyInFather( Box_1, Face_6, 'Face_6' )
geompy.addToStudyInFather( Face_2, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_2, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_2, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Face_2, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Face_3, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_3, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_3, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Face_3, Edge_4, 'Edge_4' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Box_1)
Prism_3D = Mesh_1.Prism()
Regular_1D = Mesh_1.Segment(geom=Edge_2)
Number_of_Segments_1 = Regular_1D.NumberOfSegments("mesh_nx",None,[])
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Regular_1D_1 = Mesh_1.Segment(geom=Edge_1)
Number_of_Segments_2 = Regular_1D_1.NumberOfSegments("mesh_ny",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_1)
Regular_1D_2 = Mesh_1.Segment(geom=Edge_5)
Number_of_Segments_3 = Regular_1D_2.NumberOfSegments("mesh_nz",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_5)
Regular_1D_3 = Mesh_1.Segment()
isDone = Mesh_1.Compute()
Group_1_0 = Mesh_1.GroupOnGeom(Face_1,'Group_1_0',SMESH.FACE)
[ Group_1_0 ] = Mesh_1.GetGroups()
Group_2_0 = Mesh_1.GroupOnGeom(Face_6,'Group_2_0',SMESH.FACE)
[ Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Group_3_0 = Mesh_1.GroupOnGeom(Face_2,'Group_3_0',SMESH.FACE)
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
Group_4_0 = Mesh_1.GroupOnGeom(Face_5,'Group_4_0',SMESH.FACE)
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0 ] = Mesh_1.GetGroups()
Group_5_0 = Mesh_1.GroupOnGeom(Face_3,'Group_5_0',SMESH.FACE)
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0, Group_5_0 ] = Mesh_1.GetGroups()
Group_6_0 = Mesh_1.GroupOnGeom(Face_4,'Group_6_0',SMESH.FACE)
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0, Group_5_0, Group_6_0 ] = Mesh_1.GetGroups()
Sub_mesh_1 = Regular_1D.GetSubMesh()
Sub_mesh_2 = Regular_1D_1.GetSubMesh()
Sub_mesh_3 = Regular_1D_2.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Prism_3D.GetAlgorithm(), 'Prism_3D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_5_0, 'Group_5_0')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Group_6_0, 'Group_6_0')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_2')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
