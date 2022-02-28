#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
<<<<<<< HEAD
sys.path.insert(0, r'/home/student/software/femus/applications/2021_fall/ayillah/input')
=======
sys.path.insert(0, r'/home/gbornia/software/femus/applications/2021_fall/gbornia/input')
>>>>>>> femttu/master

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("l_x", 3)
notebook.set("l_y", 1)
notebook.set("l_x_half", "0.5*l_x")
notebook.set("l_y_half", "0.5*l_y")
<<<<<<< HEAD
notebook.set("l_x", 3)
notebook.set("l_y", 1)
notebook.set("l_x_half", "0.5*l_x")
notebook.set("l_y_half", "0.5*l_y")
notebook.set("n_x", 4)
notebook.set("n_y", 3)
=======
>>>>>>> femttu/master
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
Face_1 = geompy.MakeFaceHW("l_x", "l_y", 1)
Translation_1 = geompy.MakeTranslation(Face_1, "l_x_half", "l_y_half", 0)
<<<<<<< HEAD
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Translation_1, geompy.ShapeType["EDGE"], True)
[Edge_1, Edge_2, Edge_3, Edge_4] = geompy.GetExistingSubObjects(Translation_1, False)
=======
>>>>>>> femttu/master
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
<<<<<<< HEAD
geompy.addToStudyInFather( Translation_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Translation_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Translation_1, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Translation_1, Edge_4, 'Edge_4' )
=======
>>>>>>> femttu/master

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Translation_1)
<<<<<<< HEAD
Edge_1_1 = Mesh_1.GroupOnGeom(Edge_1,'Edge_1',SMESH.EDGE)
Edge_2_1 = Mesh_1.GroupOnGeom(Edge_2,'Edge_2',SMESH.EDGE)
Edge_3_1 = Mesh_1.GroupOnGeom(Edge_3,'Edge_3',SMESH.EDGE)
Edge_4_1 = Mesh_1.GroupOnGeom(Edge_4,'Edge_4',SMESH.EDGE)
Regular_1D = Mesh_1.Segment(geom=Edge_2)
Number_of_Segments_1 = Regular_1D.NumberOfSegments("n_x")
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Regular_1D_1 = Mesh_1.Segment(geom=Edge_1)
Number_of_Segments_2 = Regular_1D_1.NumberOfSegments("n_y")
Propagation_of_1D_Hyp_1 = Regular_1D_1.Propagation()
Regular_1D_2 = Mesh_1.Segment()
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_1.Compute()
[ Edge_1_1, Edge_2_1, Edge_3_1, Edge_4_1 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
[ Edge_1_1, Edge_2_1, Edge_3_1, Edge_4_1 ] = Mesh_1.GetGroups()
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_1_0' )
nbAdd = Group_1_0.Add( [ 4, 5, 6, 7 ] )
[ Edge_1_1, Edge_2_1, Edge_3_1, Edge_4_1, Group_1_0 ] = Mesh_1.GetGroups()
Group_2_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_2_0' )
nbAdd = Group_2_0.Add( [ 1, 2, 3 ] )
[ Edge_1_1, Edge_2_1, Edge_3_1, Edge_4_1, Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Group_3_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_3_0' )
nbAdd = Group_3_0.Add( [ 11, 12, 13, 14 ] )
[ Edge_1_1, Edge_2_1, Edge_3_1, Edge_4_1, Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
Group_4_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_4_0' )
nbAdd = Group_4_0.Add( [ 8, 9, 10 ] )
[ Edge_1_1, Edge_2_1, Edge_3_1, Edge_4_1, Group_1_0, Group_2_0, Group_3_0, Group_4_0 ] = Mesh_1.GetGroups()
Sub_mesh_1 = Regular_1D.GetSubMesh()
Sub_mesh_2 = Regular_1D_1.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Edge_1_1, 'Edge_1')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Propagation_of_1D_Hyp_1, 'Propagation of 1D Hyp. on Opposite Edges_2')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Edge_4_1, 'Edge_4')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Edge_2_1, 'Edge_2')
smesh.SetName(Edge_3_1, 'Edge_3')
=======
NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_2D)
isDone = Mesh_1.Compute()


## Set names of Mesh objects
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
>>>>>>> femttu/master


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
