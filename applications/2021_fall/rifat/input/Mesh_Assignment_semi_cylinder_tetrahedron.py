#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/student/software/femus/applications/2021_fall/rifat/input')

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
Cylinder_1 = geompy.MakeCylinderRH(1, 2)
Box_1 = geompy.MakeBoxDXDYDZ(2, 2, 2)
Translation_1 = geompy.MakeTranslation(Box_1, -1, 0, 0)
Cut_1 = geompy.MakeCutList(Cylinder_1, [Translation_1], True)
[Face_1,Face_2,Face_3,Face_4] = geompy.ExtractShapes(Cut_1, geompy.ShapeType["FACE"], True)
[Face_1, Face_2, Face_3, Face_4] = geompy.GetExistingSubObjects(Cut_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudyInFather( Cut_1, Face_1, 'Face_1' )
geompy.addToStudyInFather( Cut_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Cut_1, Face_3, 'Face_3' )
geompy.addToStudyInFather( Cut_1, Face_4, 'Face_4' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Cut_1)
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 0.6 )
NETGEN_3D_Parameters_1.SetMinSize( 0.25 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 2 )
NETGEN_3D_Parameters_1.SetChordalError( -1 )
NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_1.SetCheckChartBoundary( 80 )
Viscous_Layers_1 = NETGEN_1D_2D_3D.ViscousLayers(1,1,1,[],1,smeshBuilder.SURF_OFFSET_SMOOTH)
isDone = Mesh_1.Compute()
Group_1_0 = Mesh_1.GroupOnGeom(Face_1,'Face_1',SMESH.FACE)
Group_2_0 = Mesh_1.GroupOnGeom(Face_2,'Face_2',SMESH.FACE)
Group_3_0 = Mesh_1.GroupOnGeom(Face_3,'Face_3',SMESH.FACE)
Group_4_0 = Mesh_1.GroupOnGeom(Face_4,'Face_4',SMESH.FACE)
Group_1_0.SetName( 'Group_1_0' )
Group_2_0.SetName( 'Group_2_0' )
Group_3_0.SetName( 'Group_3_0' )
Group_4_0.SetName( 'Group_4_0' )
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/rifat/input/Mesh_Assignment_semi_cylinder_tetrahedron.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Viscous_Layers_1, 'Viscous Layers_1')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_4_0, 'Group_4_0')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
