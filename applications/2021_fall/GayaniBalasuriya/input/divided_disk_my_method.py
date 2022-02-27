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
nsubs_radial = Regular_1D_1.NumberOfSegments(3,None,[])
Regular_1D_2 = Mesh_1.Segment(geom=Edge_3)
nsubs_central_ymx = Regular_1D_2.NumberOfSegments(4,None,[])
#Mesh_1.GetMesh().RemoveSubMesh( smeshObj_1 ) ### smeshObj_1 has not been yet created
#Mesh_1.GetMesh().RemoveSubMesh( smeshObj_2 ) ### smeshObj_2 has not been yet created
placeholder = Regular_1D.NumberOfSegments(15)
Regular_1D_3 = Mesh_1.Segment(geom=Edge_1)
status = Mesh_1.AddHypothesis(nsubs_radial,Edge_1)
Propagation_of_1D_Hyp = Regular_1D_1.Propagation()
Regular_1D_4 = Mesh_1.Segment(geom=Edge_2)
status = Mesh_1.AddHypothesis(nsubs_central_ymx,Edge_2)
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_2)
isDone = Mesh_1.Compute()
Mesh_2 = smesh.Mesh(Face_2)
status = Mesh_2.AddHypothesis(placeholder)
Regular_1D_5 = Mesh_2.Segment()
Quadrangle_2D_1 = Mesh_2.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_6 = Mesh_2.Segment(geom=Edge_5)
status = Mesh_2.AddHypothesis(nsubs_radial,Edge_5)
status = Mesh_2.AddHypothesis(Propagation_of_1D_Hyp,Edge_5)
Regular_1D_7 = Mesh_2.Segment(geom=Edge_6)
nsubs_central_ypx = Regular_1D_7.NumberOfSegments(5)
status = Mesh_2.AddHypothesis(Propagation_of_1D_Hyp,Edge_6)
isDone = Mesh_2.Compute()
Mesh_3 = smesh.Mesh(Face_5)
status = Mesh_3.AddHypothesis(placeholder)
Regular_1D_8 = Mesh_3.Segment()
Quadrangle_2D_2 = Mesh_3.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_9 = Mesh_3.Segment(geom=Edge_17)
status = Mesh_3.AddHypothesis(nsubs_radial,Edge_17)
status = Mesh_3.AddHypothesis(Propagation_of_1D_Hyp,Edge_17)
Regular_1D_10 = Mesh_3.Segment(geom=Edge_19)
status = Mesh_3.AddHypothesis(nsubs_central_ymx,Edge_19)
status = Mesh_3.AddHypothesis(Propagation_of_1D_Hyp,Edge_19)
isDone = Mesh_3.Compute()
Mesh_4 = smesh.Mesh(Face_4)
status = Mesh_4.AddHypothesis(placeholder)
Regular_1D_11 = Mesh_4.Segment()
Quadrangle_2D_3 = Mesh_4.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_12 = Mesh_4.Segment(geom=Edge_13)
status = Mesh_4.AddHypothesis(nsubs_radial,Edge_13)
status = Mesh_4.AddHypothesis(Propagation_of_1D_Hyp,Edge_13)
Regular_1D_13 = Mesh_4.Segment(geom=Edge_15)
status = Mesh_4.AddHypothesis(nsubs_central_ypx,Edge_15)
status = Mesh_4.AddHypothesis(Propagation_of_1D_Hyp,Edge_15)
isDone = Mesh_4.Compute()
Mesh_5 = smesh.Mesh(Face_3)
status = Mesh_5.AddHypothesis(placeholder)
Regular_1D_14 = Mesh_5.Segment()
Quadrangle_2D_4 = Mesh_5.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_15 = Mesh_5.Segment(geom=Edge_9)
status = Mesh_5.AddHypothesis(nsubs_central_ymx,Edge_9)
Regular_1D_16 = Mesh_5.Segment(geom=Edge_10)
status = Mesh_5.AddHypothesis(nsubs_central_ypx,Edge_10)
status = Mesh_5.AddHypothesis(Propagation_of_1D_Hyp,Edge_10)
status = Mesh_5.AddHypothesis(Propagation_of_1D_Hyp,Edge_9)
isDone = Mesh_5.Compute()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Mesh_2.ConvertToQuadratic(0, Mesh_2,True)
Mesh_3.ConvertToQuadratic(0, Mesh_3,True)
Mesh_4.ConvertToQuadratic(0, Mesh_4,True)
Mesh_5.ConvertToQuadratic(0, Mesh_5,True)
Mesh_compound = smesh.Concatenate( [ Mesh_1.GetMesh(), Mesh_2.GetMesh(), Mesh_3.GetMesh(), Mesh_4.GetMesh(), Mesh_5.GetMesh() ], 1, 1, 1e-05, False )
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.EDGE,SMESH.FT_FreeBorders,SMESH.FT_Undefined,0,SMESH.FT_LogicalNOT)
aCriteria.append(aCriterion)
isDone = Mesh_compound.RemoveElements( [ 1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 27, 28, 29, 35, 36, 37, 38, 39, 55, 56, 57, 62, 63, 64, 65, 83, 84, 85, 86, 87 ] )
Group_1_0 = Mesh_compound.CreateEmptyGroup( SMESH.EDGE, 'Group_1_0' )
nbAdd = Group_1_0.AddFrom( Mesh_compound.GetMesh() )
smesh.SetName(Mesh_compound, 'Mesh_compound')
try:
  Mesh_compound.ExportMED(r'/home/gbornia/software/femus/applications/OptimalControl/00_laplacian_all_dims/input/disk_quad.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Sub_mesh_1 = Regular_1D_3.GetSubMesh()
Sub_mesh_2 = Regular_1D_4.GetSubMesh()
Sub_mesh_3 = Regular_1D_6.GetSubMesh()
Sub_mesh_4 = Regular_1D_7.GetSubMesh()
Sub_mesh_5 = Regular_1D_9.GetSubMesh()
Sub_mesh_6 = Regular_1D_10.GetSubMesh()
Sub_mesh_7 = Regular_1D_12.GetSubMesh()
Sub_mesh_8 = Regular_1D_13.GetSubMesh()
Sub_mesh_9 = Regular_1D_15.GetSubMesh()
Sub_mesh_10 = Regular_1D_16.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(nsubs_central_ymx, 'nsubs_central_ymx')
smesh.SetName(Sub_mesh_8, 'Sub-mesh_8')
smesh.SetName(nsubs_radial, 'nsubs_radial')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Sub_mesh_7, 'Sub-mesh_7')
smesh.SetName(nsubs_central_ypx, 'nsubs_central_ypx')
smesh.SetName(placeholder, 'placeholder')
smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
smesh.SetName(Mesh_compound.GetMesh(), 'Mesh_compound')
smesh.SetName(Sub_mesh_6, 'Sub-mesh_6')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Mesh_3.GetMesh(), 'Mesh_3')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')
smesh.SetName(Mesh_5.GetMesh(), 'Mesh_5')
smesh.SetName(Mesh_4.GetMesh(), 'Mesh_4')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Sub_mesh_10, 'Sub-mesh_10')
smesh.SetName(Sub_mesh_9, 'Sub-mesh_9')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
