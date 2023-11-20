#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.8.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("lx", 1)
notebook.set("ly", 1)
notebook.set("origin_x", 0)
notebook.set("origin_y", 0)
notebook.set("vertex_x", "origin_x+ lx")
notebook.set("vertex_y", "origin_y + ly")
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
Vertex_1 = geompy.MakeVertex("origin_x", "origin_y", 0)
Vertex_2 = geompy.MakeVertex("vertex_x", "origin_y", 0)
Vertex_3 = geompy.MakeVertex("origin_x", "vertex_y", 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_3 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_1)
Face_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_3], 1)
[Edge_1,Edge_2,Edge_3] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)
[Edge_1, Edge_2, Edge_3] = geompy.GetExistingSubObjects(Face_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudyInFather( Face_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_1, Edge_3, 'Edge_3' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

aFilterManager = smesh.CreateFilterManager()
Mesh_1 = smesh.Mesh(Face_1)
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Simple_Parameters_1 = NETGEN_1D_2D.Parameters(smeshBuilder.SIMPLE)
NETGEN_2D_Simple_Parameters_1.SetNumberOfSegments( 1 )
NETGEN_2D_Simple_Parameters_1.LengthFromEdges()
NETGEN_2D_Simple_Parameters_1.SetAllowQuadrangles( 0 )
Group_1_0 = Mesh_1.GroupOnGeom(Edge_1,'Edge_1',SMESH.EDGE)
Group_2_0 = Mesh_1.GroupOnGeom(Edge_2,'Edge_2',SMESH.EDGE)
Group_3_0 = Mesh_1.GroupOnGeom(Edge_3,'Edge_3',SMESH.EDGE)
isDone = Mesh_1.Compute()
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0)
Group_1_0.SetName( 'Group_1_0' )
Group_2_0.SetName( 'Group_2_0' )
Group_3_0.SetName( 'Group_3_0' )
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,Face_1)
aCriteria.append(aCriterion)
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,'Face_1')
aCriteria.append(aCriterion)
aFilter_2 = smesh.GetFilterFromCriteria(aCriteria)
aFilter_2.SetMesh(Mesh_1.GetMesh())
aFilter0x4fd71d0 = aFilterManager.CreateFilter()
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,'Face_1')
aCriteria.append(aCriterion)
aFilter_2.SetCriteria(aCriteria)
Group_4_0 = Mesh_1.GroupOnFilter( SMESH.FACE, 'Group_1', aFilter_2 )
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0 ] = Mesh_1.GetGroups()
Group_4_0.SetName( 'Group_4_0' )
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')

Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')

## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(NETGEN_2D_Simple_Parameters_1, 'NETGEN 2D Simple Parameters_1')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_2_0, 'Group_2_0')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
