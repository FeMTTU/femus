#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/applications/2021_fall/giorgiobornia/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("l_x", 1)
notebook.set("l_y", 1)
notebook.set("mesh_n_x", 2)
notebook.set("mesh_n_y", 3)
notebook.set("half_l_x", "0.5*l_x")
notebook.set("half_l_y", "0.5*l_y")
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
[geomObj_1,geomObj_2,geomObj_3,geomObj_4] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)
Translation_1 = geompy.MakeTranslation(Face_1, "half_l_x", "half_l_y", 0)
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Translation_1, geompy.ShapeType["EDGE"], True)
[Edge_5, Edge_6, Edge_7, Edge_8] = geompy.GetExistingSubObjects(Translation_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudyInFather( Translation_1, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Translation_1, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Translation_1, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Translation_1, Edge_8, 'Edge_8' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Translation_1)
Regular_1D = Mesh_1.Segment()
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Edge_5_1 = Mesh_1.GroupOnGeom(Edge_5,'Edge_5',SMESH.EDGE)
Edge_6_1 = Mesh_1.GroupOnGeom(Edge_6,'Edge_6',SMESH.EDGE)
Edge_7_1 = Mesh_1.GroupOnGeom(Edge_7,'Edge_7',SMESH.EDGE)
Edge_8_1 = Mesh_1.GroupOnGeom(Edge_8,'Edge_8',SMESH.EDGE)
Number_of_Segments_1 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_1.SetNumberOfSegments( "mesh_n_x" )
Regular_1D_1 = Mesh_1.Segment(geom=Edge_5)
Number_of_Segments_2 = Regular_1D_1.NumberOfSegments("mesh_n_y",None,[])
Propagation_of_1D_Hyp = Regular_1D_1.Propagation()
[ Edge_5_1, Edge_6_1, Edge_7_1, Edge_8_1 ] = Mesh_1.GetGroups()
Regular_1D_2 = Mesh_1.Segment(geom=Edge_6)
Number_of_Segments_3 = Regular_1D_2.NumberOfSegments("mesh_n_x",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_6)
isDone = Mesh_1.Compute()
[ Edge_5_1, Edge_6_1, Edge_7_1, Edge_8_1 ] = Mesh_1.GetGroups()
Mesh_1.QuadTo4Tri( Mesh_1 )
[ Edge_5_1, Edge_6_1, Edge_7_1, Edge_8_1 ] = Mesh_1.GetGroups()
Sub_mesh_1 = Regular_1D_1.GetSubMesh()
Sub_mesh_2 = Regular_1D_2.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Edge_5_1, 'Edge_5')
smesh.SetName(Edge_6_1, 'Edge_6')
smesh.SetName(Edge_7_1, 'Edge_7')
smesh.SetName(Edge_8_1, 'Edge_8')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_2')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
