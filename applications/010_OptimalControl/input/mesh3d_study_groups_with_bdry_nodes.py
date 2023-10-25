#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.8.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/applications/OptimalControl/01_boundary_control_inequality/00_dirichlet/00_dirichlet_boundary_fractional/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("lx", 1)
notebook.set("ly", 1)
notebook.set("lz", 1)
notebook.set("nx", 2)
notebook.set("ny", 4)
notebook.set("nz", 1)
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
Box_1 = geompy.MakeBoxDXDYDZ("lx", "ly", "lz")
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
Regular_1D = Mesh_1.Segment(geom=Edge_2)
Regular_1D_1 = Mesh_1.Segment(geom=Edge_3)
Regular_1D_2 = Mesh_1.Segment(geom=Edge_5)
Cartesian_3D = smesh.CreateHypothesis('Cartesian_3D')
Regular_1D_3 = Mesh_1.Segment(geom=Edge_1)
Regular_1D_4 = Mesh_1.Segment()
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Hexa_3D = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
Number_of_Segments_7 = Regular_1D.NumberOfSegments("nx",None,[])
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Number_of_Segments_8 = Regular_1D_2.NumberOfSegments("nz",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_5)
Number_of_Segments_9 = Regular_1D_3.NumberOfSegments("ny",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_1)
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Group_1_0 = Mesh_1.GroupOnGeom(Face_1,'Group_1_0',SMESH.FACE)
Group_2_0 = Mesh_1.GroupOnGeom(Face_6,'Group_2_0',SMESH.FACE)
Group_3_0 = Mesh_1.GroupOnGeom(Face_2,'Group_3_0',SMESH.FACE)
Group_4_0 = Mesh_1.GroupOnGeom(Face_5,'Group_4_0',SMESH.FACE)
Group_5_0 = Mesh_1.GroupOnGeom(Face_3,'Group_5_0',SMESH.FACE)
Group_6_0 = Mesh_1.GroupOnGeom(Face_4,'Group_6_0',SMESH.FACE)
Group_7_0 = Mesh_1.CreateEmptyGroup( SMESH.NODE, 'Group_7_0' )
nbAdd = Group_7_0.Add( [ 18, 48, 19, 49, 20, 70, 17, 44, 16, 43, 15, 66 ] )
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/gbornia/software/femus_build/applications/OptimalControl/01_boundary_control_inequality/00_dirichlet/00_dirichlet_boundary_fractional/input/Mesh_3_groups_with_bdry_nodes_coarser.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Sub_mesh_1 = Regular_1D.GetSubMesh()
Regular_1D_5 = Regular_1D_1.GetSubMesh()
Sub_mesh_3 = Regular_1D_2.GetSubMesh()
Sub_mesh_2 = Regular_1D_3.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Group_7_0, 'Group_7_0')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Cartesian_3D, 'Cartesian_3D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_5_0, 'Group_5_0')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Group_6_0, 'Group_6_0')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Number_of_Segments_8, 'Number of Segments_8')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_7')
smesh.SetName(Number_of_Segments_7, 'Number of Segments_7')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Number_of_Segments_9, 'Number of Segments_9')
smesh.SetName(Regular_1D_5, 'Regular_1D')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
