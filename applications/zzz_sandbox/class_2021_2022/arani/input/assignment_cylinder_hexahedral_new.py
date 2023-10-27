#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/student/software/femus/applications/2021_fall/arani/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("geom_r", 1)
notebook.set("geom_h", 1)
notebook.set("mesh_n_theta_1", 2)
notebook.set("mesh_n_theta_2", 2)
notebook.set("mesh_n_radial", 2)
notebook.set("mesh_n_h", 2)
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
Divided_Cylinder_1 = geompy.MakeDividedCylinder("geom_r", "geom_h", GEOM.SQUARE)
Translation_1 = geompy.MakeTranslation(Divided_Cylinder_1, "geom_r", "geom_r", 0)
[Face_1,Face_2,Face_3,Face_4,Face_5,Face_6,Face_7,Face_8,Face_9,Face_10,Face_11,Face_12,Face_13,Face_14,Face_15,Face_16,Face_17,Face_18,Face_19,Face_20,Face_21,Face_22] = geompy.ExtractShapes(Translation_1, geompy.ShapeType["FACE"], True)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_2, geompy.ShapeType["EDGE"], True)
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Face_4, geompy.ShapeType["EDGE"], True)
[Edge_9,Edge_10,Edge_11,Edge_12] = geompy.ExtractShapes(Face_16, geompy.ShapeType["EDGE"], True)
[Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, Face_7, Face_8, Face_9, Face_10, Face_11, Face_12, Face_13, Face_14, Face_15, Face_16, Face_17, Face_18, Face_19, Face_20, Face_21, Face_22] = geompy.GetExistingSubObjects(Translation_1, False)

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Translation_1)
Prism_3D = Mesh_1.Prism()
Lower_face1 = Mesh_1.GroupOnGeom(Face_4,'Face_4',SMESH.FACE)
Face_5_1 = Mesh_1.GroupOnGeom(Face_5,'Face_5',SMESH.FACE)
Lower_face2 = Mesh_1.GroupOnGeom(Face_6,'Face_6',SMESH.FACE)
Face_7_1 = Mesh_1.GroupOnGeom(Face_7,'Face_7',SMESH.FACE)
Lower_face3 = Mesh_1.GroupOnGeom(Face_11,'Face_11',SMESH.FACE)
Face_12_1 = Mesh_1.GroupOnGeom(Face_12,'Face_12',SMESH.FACE)
Lower_face4 = Mesh_1.GroupOnGeom(Face_16,'Face_16',SMESH.FACE)
Face_17_1 = Mesh_1.GroupOnGeom(Face_17,'Face_17',SMESH.FACE)
Lower_face5 = Mesh_1.GroupOnGeom(Face_18,'Face_18',SMESH.FACE)
Cyl_face3 = Mesh_1.GroupOnGeom(Face_20,'Face_20',SMESH.FACE)
Cyl_face4 = Mesh_1.GroupOnGeom(Face_21,'Face_21',SMESH.FACE)
Group_3 = Mesh_1.GroupOnGeom(Face_22,'Face_22',SMESH.FACE)
Regular_1D = Mesh_1.Segment(geom=Edge_1)
Number_of_Segments_1 = Regular_1D.NumberOfSegments("mesh_n_h",None,[])
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Number_of_Segments_2 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_2.SetNumberOfSegments( "mesh_n_radial" )
Regular_1D_1 = Mesh_1.Segment(geom=Edge_5)
Number_of_Segments_3 = Regular_1D_1.NumberOfSegments("mesh_n_radial",None,[])
Regular_1D_2 = Mesh_1.Segment(geom=Edge_2)
Number_of_Segments_4 = Regular_1D_2.NumberOfSegments("mesh_n_theta_1",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_2)
Regular_1D_3 = Mesh_1.Segment(geom=Edge_11)
Number_of_Segments_5 = Regular_1D_3.NumberOfSegments("mesh_n_theta_2",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_11)
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_5)
[ smeshObj_1, smeshObj_2, smeshObj_3, Lower_face1, Face_5_1, Lower_face2, Face_7_1, smeshObj_4, smeshObj_5, smeshObj_6, Lower_face3, Face_12_1, smeshObj_7, smeshObj_8, smeshObj_9, Lower_face4, Face_17_1, Lower_face5, smeshObj_10, Cyl_face3, Cyl_face4, Group_3 ] = Mesh_1.GetGroups()
Regular_1D_4 = Mesh_1.Segment()
[ smeshObj_1, smeshObj_2, smeshObj_3, Lower_face1, Face_5_1, Lower_face2, Face_7_1, smeshObj_4, smeshObj_5, smeshObj_6, Lower_face3, Face_12_1, smeshObj_7, smeshObj_8, smeshObj_9, Lower_face4, Face_17_1, Lower_face5, smeshObj_10, Cyl_face3, Cyl_face4, Group_3 ] = Mesh_1.GetGroups()
isDone = Mesh_1.Compute()
[ smeshObj_1, smeshObj_2, smeshObj_3, Lower_face1, Face_5_1, Lower_face2, Face_7_1, smeshObj_4, smeshObj_5, smeshObj_6, Lower_face3, Face_12_1, smeshObj_7, smeshObj_8, smeshObj_9, Lower_face4, Face_17_1, Lower_face5, smeshObj_10, Cyl_face3, Cyl_face4, Group_3 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
[ smeshObj_2, smeshObj_3, Lower_face1, Face_5_1, Lower_face2, Face_7_1, Lower_face3, Face_12_1, Lower_face4, Face_17_1, Lower_face5, smeshObj_10, Cyl_face3, Cyl_face4 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/arani/input/assignment_mesh_cylinder_hexahedral_new.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Cyl_face3.SetName( 'Cyl_face3' )
Cyl_face4.SetName( 'Cyl_face4' )
Group_3 = Mesh_1.GetMesh().UnionListOfGroups([ Cyl_face3, Cyl_face4 ], 'Group_3' )
Group_3_2 = Mesh_1.GetMesh().UnionListOfGroups([ smeshObj_2, smeshObj_3 ], 'Group_3_2' )
Group_3_0 = Mesh_1.GetMesh().UnionListOfGroups([ Group_3, Group_3_2 ], 'Group_3_0' )
Lower_face1.SetName( 'Lower_face1' )
Lower_face2.SetName( 'Lower_face2' )
Lower_face3.SetName( 'Lower_face_3' )
Lower_face4.SetName( 'Lower_face4' )
Lower_face5.SetName( 'Lower_face5' )
Lower_face3.SetName( 'Lower_face3' )
Group_1_0 = Mesh_1.GetMesh().UnionListOfGroups([ Lower_face1, Lower_face2, Lower_face3, Lower_face4, Lower_face5 ], 'Group_1_0' )
Group_2_0 = Mesh_1.GetMesh().UnionListOfGroups([ Face_5_1, Face_7_1, Face_12_1, Face_17_1, smeshObj_10 ], 'Group_2_0' )
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/arani/input/assignment_mesh_cylinder_hexahedral_new.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Sub_mesh_1 = Regular_1D.GetSubMesh()
Sub_mesh_2 = Regular_1D_1.GetSubMesh()
Sub_mesh_4 = Regular_1D_3.GetSubMesh()
Sub_mesh_3 = Mesh_1.GetSubMesh( Edge_6, 'Sub-mesh_3' )

## some objects were removed
aStudyBuilder = salome.myStudy.NewBuilder()
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_10))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_6))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_7))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_8))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_9))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_5))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_4))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_1))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_3))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_2))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)

## Set names of Mesh objects
smesh.SetName(Lower_face5, 'Lower_face5')
smesh.SetName(Lower_face3, 'Lower_face3')
smesh.SetName(Face_12_1, 'Face_12')
smesh.SetName(Lower_face4, 'Lower_face4')
smesh.SetName(Face_17_1, 'Face_17')
smesh.SetName(Prism_3D.GetAlgorithm(), 'Prism_3D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Face_5_1, 'Face_5')
smesh.SetName(Lower_face1, 'Lower_face1')
smesh.SetName(Face_7_1, 'Face_7')
smesh.SetName(Lower_face2, 'Lower_face2')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_3_2, 'Group_3_2')
smesh.SetName(Group_3, 'Group_3')
smesh.SetName(Cyl_face4, 'Cyl_face4')
smesh.SetName(Cyl_face3, 'Cyl_face3')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_4')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_5, 'Number of Segments_5')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
