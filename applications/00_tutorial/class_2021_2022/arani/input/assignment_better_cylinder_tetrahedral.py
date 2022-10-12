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
notebook.set("r", 1)
notebook.set("z", 2)
notebook.set("n_seg", 4)
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
Cylinder_1 = geompy.MakeCylinderRH("r", "z")
Translation_1 = geompy.MakeTranslation(Cylinder_1, "r", "r", 0)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Translation_1)
NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_2D)
NETGEN_3D = Mesh_1.Tetrahedron()
CompositeSegment_1D = Mesh_1.Segment(algo=smeshBuilder.COMPOSITE)
Number_of_Segments_1 = CompositeSegment_1D.NumberOfSegments("n_seg",None,[])
status = Mesh_1.RemoveHypothesis(CompositeSegment_1D)
Regular_1D = Mesh_1.Segment()
Number_of_Segments_1.SetNumberOfSegments( "n_seg" )
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/arani/input/assignment_mesh_cylinder_tetrahedral3.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
NETGEN_2D_Parameters_1 = NETGEN_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( 1 )
NETGEN_2D_Parameters_1.SetMinSize( 0.1 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 2 )
NETGEN_2D_Parameters_1.SetChordalError( -1 )
NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetWorstElemMeasure( 0 )
NETGEN_2D_Parameters_1.SetUseDelauney( 0 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_1.SetCheckChartBoundary( 32 )
NETGEN_3D_Parameters_1 = NETGEN_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 1 )
NETGEN_3D_Parameters_1.SetMinSize( 0.1 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 2 )
NETGEN_3D_Parameters_1.SetElemSizeWeight( 3.98838e-43 )
NETGEN_3D_Parameters_1.SetNbVolOptSteps( 4 )
NETGEN_3D_Parameters_1.SetCheckOverlapping( 0 )
NETGEN_3D_Parameters_1.SetCheckChartBoundary( 128 )
NETGEN_3D_Parameters_1_1 = smesh.CreateHypothesis('NETGEN_Parameters', 'NETGENEngine')
NETGEN_3D_Parameters_1_1.SetMaxSize( 10.3464 )
NETGEN_3D_Parameters_1_1.SetMinSize( 0 )
NETGEN_3D_Parameters_1_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1_1.SetFineness( 2 )
NETGEN_3D_Parameters_1_1.SetChordalError( -1 )
NETGEN_3D_Parameters_1_1.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_1_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1_1.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_1_1.SetCheckChartBoundary( 208 )
status = Mesh_1.RemoveHypothesis(Regular_1D)
status = Mesh_1.RemoveHypothesis(NETGEN_2D)
status = Mesh_1.RemoveHypothesis(NETGEN_3D)
status = Mesh_1.RemoveHypothesis(Number_of_Segments_1)
status = Mesh_1.RemoveHypothesis(NETGEN_2D_Parameters_1)
status = Mesh_1.RemoveHypothesis(NETGEN_3D_Parameters_1)
NETGEN_3D_Parameters_2 = smesh.CreateHypothesis('NETGEN_Parameters', 'NETGENEngine')
NETGEN_3D_Parameters_2.SetMaxSize( 1 )
NETGEN_3D_Parameters_2.SetMinSize( 0.1 )
NETGEN_3D_Parameters_2.SetSecondOrder( 0 )
NETGEN_3D_Parameters_2.SetOptimize( 1 )
NETGEN_3D_Parameters_2.SetFineness( 2 )
NETGEN_3D_Parameters_2.SetChordalError( -1 )
NETGEN_3D_Parameters_2.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_2.SetFuseEdges( 1 )
NETGEN_3D_Parameters_2.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_2.SetCheckChartBoundary( 224 )
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Simple_Parameters_1 = NETGEN_1D_2D_3D.Parameters(smeshBuilder.SIMPLE)
NETGEN_3D_Simple_Parameters_1.SetNumberOfSegments( 6 )
NETGEN_3D_Simple_Parameters_1.LengthFromEdges()
NETGEN_3D_Simple_Parameters_1.LengthFromFaces()
NETGEN_3D_Parameters_3 = smesh.CreateHypothesis('NETGEN_Parameters', 'NETGENEngine')
NETGEN_3D_Parameters_3.SetMaxSize( 1 )
NETGEN_3D_Parameters_3.SetMinSize( 0 )
NETGEN_3D_Parameters_3.SetSecondOrder( 0 )
NETGEN_3D_Parameters_3.SetOptimize( 1 )
NETGEN_3D_Parameters_3.SetFineness( 2 )
NETGEN_3D_Parameters_3.SetChordalError( -1 )
NETGEN_3D_Parameters_3.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_3.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_3.SetFuseEdges( 1 )
NETGEN_3D_Parameters_3.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_3.SetCheckChartBoundary( 96 )
status = Mesh_1.RemoveHypothesis(NETGEN_3D_Simple_Parameters_1)
NETGEN_3D_Simple_Parameters_2 = NETGEN_1D_2D_3D.Parameters(smeshBuilder.SIMPLE)
NETGEN_3D_Simple_Parameters_2.SetNumberOfSegments( "n_seg" )
NETGEN_3D_Simple_Parameters_2.LengthFromEdges()
NETGEN_3D_Simple_Parameters_2.LengthFromFaces()
NETGEN_3D_Parameters_4 = smesh.CreateHypothesis('NETGEN_Parameters', 'NETGENEngine')
NETGEN_3D_Parameters_4.SetMaxSize( 0.6 )
NETGEN_3D_Parameters_4.SetMinSize( 0.5 )
NETGEN_3D_Parameters_4.SetSecondOrder( 0 )
NETGEN_3D_Parameters_4.SetOptimize( 1 )
NETGEN_3D_Parameters_4.SetFineness( 2 )
NETGEN_3D_Parameters_4.SetChordalError( -1 )
NETGEN_3D_Parameters_4.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_4.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_4.SetFuseEdges( 1 )
NETGEN_3D_Parameters_4.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_4.SetCheckChartBoundary( 160 )
status = Mesh_1.RemoveHypothesis(NETGEN_3D_Simple_Parameters_2)
NETGEN_3D_Simple_Parameters_3 = NETGEN_1D_2D_3D.Parameters(smeshBuilder.SIMPLE)
NETGEN_3D_Simple_Parameters_3.SetNumberOfSegments( "n_seg" )
NETGEN_3D_Simple_Parameters_3.LengthFromEdges()
NETGEN_3D_Simple_Parameters_3.LengthFromFaces()
status = Mesh_1.RemoveHypothesis(NETGEN_3D_Simple_Parameters_3)
NETGEN_3D_Simple_Parameters_4 = NETGEN_1D_2D_3D.Parameters(smeshBuilder.SIMPLE)
NETGEN_3D_Simple_Parameters_4.SetNumberOfSegments( "n_seg" )
NETGEN_3D_Simple_Parameters_4.LengthFromEdges()
NETGEN_3D_Simple_Parameters_4.LengthFromFaces()
NETGEN_3D_Simple_Parameters_5 = smesh.CreateHypothesis('NETGEN_SimpleParameters_3D', 'NETGENEngine')
NETGEN_3D_Simple_Parameters_5.SetNumberOfSegments( "n_seg" )
NETGEN_3D_Simple_Parameters_5.LengthFromEdges()
NETGEN_3D_Simple_Parameters_5.LengthFromFaces()
status = Mesh_1.RemoveHypothesis(NETGEN_3D_Simple_Parameters_4)
NETGEN_3D_Parameters_5 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_5.SetMaxSize( 0.5 )
NETGEN_3D_Parameters_5.SetMinSize( 0.5 )
NETGEN_3D_Parameters_5.SetSecondOrder( 0 )
NETGEN_3D_Parameters_5.SetOptimize( 1 )
NETGEN_3D_Parameters_5.SetFineness( 2 )
NETGEN_3D_Parameters_5.SetChordalError( -1 )
NETGEN_3D_Parameters_5.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_5.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_5.SetFuseEdges( 1 )
NETGEN_3D_Parameters_5.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_5.SetCheckChartBoundary( 32 )
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0)
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/arani/input/assignment_mesh_cylinder_tetrahedral3.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/arani/input/assignment_mesh_cylinder_tetrahedral.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.FACE, 'Group_1_0' )
nbAdd = Group_1_0.Add( [ 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84 ] )
[ Group_1_0 ] = Mesh_1.GetGroups()
Group_3_0 = Mesh_1.CreateEmptyGroup( SMESH.FACE, 'Group_3_0' )
nbAdd = Group_3_0.Add( [ 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190 ] )
[ Group_1_0, Group_3_0 ] = Mesh_1.GetGroups()
Group_2_0 = Mesh_1.CreateEmptyGroup( SMESH.FACE, 'Group_2_0' )
nbAdd = Group_2_0.Add( [ 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57 ] )
[ Group_1_0, Group_3_0, Group_2_0 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/arani/input/assignment_mesh_cylinder_tetrahedral.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(NETGEN_3D_Parameters_5, 'NETGEN 3D Parameters_5')
smesh.SetName(NETGEN_3D_Simple_Parameters_4, 'NETGEN 3D Simple Parameters_4')
smesh.SetName(NETGEN_3D_Parameters_4, 'NETGEN 3D Parameters_4')
smesh.SetName(NETGEN_3D_Simple_Parameters_3, 'NETGEN 3D Simple Parameters_3')
smesh.SetName(NETGEN_3D_Simple_Parameters_5, 'NETGEN 3D Simple Parameters_5')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(NETGEN_3D.GetAlgorithm(), 'NETGEN 3D')
smesh.SetName(CompositeSegment_1D.GetAlgorithm(), 'CompositeSegment_1D')
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(NETGEN_3D_Simple_Parameters_1, 'NETGEN 3D Simple Parameters_1')
smesh.SetName(NETGEN_3D_Parameters_2, 'NETGEN 3D Parameters_2')
smesh.SetName(NETGEN_3D_Parameters_1_1, 'NETGEN 3D Parameters_1')
smesh.SetName(NETGEN_3D_Simple_Parameters_2, 'NETGEN 3D Simple Parameters_2')
smesh.SetName(NETGEN_3D_Parameters_3, 'NETGEN 3D Parameters_3')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
