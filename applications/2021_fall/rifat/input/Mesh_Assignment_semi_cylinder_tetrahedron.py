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
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/rifat/input/Mesh_Assignment_semi_cylinder_tetrahedron.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Mesh_1.ConvertToQuadratic(0)
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/rifat/input/Mesh_Assignment_semi_cylinder_tetrahedron.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Mesh_1.ConvertFromQuadratic()
Mesh_1.ConvertToQuadratic(0)
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/rifat/input/Mesh_Assignment_semi_cylinder_tetrahedron.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Mesh_2 = smesh.Mesh(Cut_1)
status = Mesh_2.AddHypothesis(Viscous_Layers_1)
NETGEN_1D_2D_3D_1 = Mesh_2.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Simple_Parameters_1 = NETGEN_1D_2D_3D_1.Parameters(smeshBuilder.SIMPLE)
NETGEN_3D_Simple_Parameters_1.SetNumberOfSegments( 10 )
NETGEN_3D_Simple_Parameters_1.LengthFromEdges()
NETGEN_3D_Simple_Parameters_1.LengthFromFaces()
NETGEN_3D_Simple_Parameters_1.SetLocalLength( 0.6 )
NETGEN_3D_Simple_Parameters_1.SetMaxElementArea( 0.3 )
NETGEN_3D_Simple_Parameters_1.SetMaxElementVolume( 0.3 )
isDone = Mesh_2.Compute()
Mesh_2.ConvertToQuadratic(0)
Group_1_0_1 = Mesh_2.CreateEmptyGroup( SMESH.FACE, 'Group_1_0' )
nbAdd = Group_1_0_1.Add( [ 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96 ] )
[ Group_1_0_1 ] = Mesh_2.GetGroups()
Group_2_0_1 = Mesh_2.CreateEmptyGroup( SMESH.FACE, 'Group_2_0' )
nbAdd = Group_2_0_1.Add( [ 23, 24, 25, 26, 27, 28, 29, 30 ] )
[ Group_1_0_1, Group_2_0_1 ] = Mesh_2.GetGroups()
Group_3_0_1 = Mesh_2.CreateEmptyGroup( SMESH.FACE, 'Group_1' )
nbAdd = Group_3_0_1.Add( [ 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58 ] )
[ Group_1_0_1, Group_2_0_1, Group_3_0_1 ] = Mesh_2.GetGroups()
Group_3_0_1.SetName( 'Group_3_0' )
[ Group_1_0_1, Group_2_0_1, Group_3_0_1 ] = Mesh_2.GetGroups()
Group_4_0_1 = Mesh_2.CreateEmptyGroup( SMESH.FACE, 'Group_1' )
nbAdd = Group_4_0_1.Add( [ 31, 32, 33, 34, 35, 36, 37, 38 ] )
[ Group_1_0_1, Group_2_0_1, Group_3_0_1, Group_4_0_1 ] = Mesh_2.GetGroups()
Group_4_0_1.SetName( 'Group_4_0' )
smesh.SetName(Mesh_2, 'Mesh_2')
try:
  Mesh_2.ExportMED(r'/home/student/software/femus/applications/2021_fall/rifat/input/Mesh_Assignment_semi_cylinder_tetrahedron.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Mesh_1.ConvertFromQuadratic()
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0 ] = Mesh_1.GetGroups()
status = Mesh_1.RemoveHypothesis(Viscous_Layers_1)
isDone = Mesh_1.Compute()
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0)
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/rifat/input/Mesh_Assignment_semi_cylinder_tetrahedron.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Viscous_Layers_1, 'Viscous Layers_1')
smesh.SetName(NETGEN_3D_Simple_Parameters_1, 'NETGEN 3D Simple Parameters_1')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Group_4_0_1, 'Group_4_0')
smesh.SetName(Group_3_0_1, 'Group_3_0')
smesh.SetName(Group_2_0_1, 'Group_2_0')
smesh.SetName(Group_1_0_1, 'Group_1_0')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
