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
notebook.set("l_x", 1)
notebook.set("l_y", 1)
notebook.set("l_z", 2)
notebook.set("l_x_p", 4)
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
Cylinder_1 = geompy.MakeCylinderRH(1, 2)
Translation_1 = geompy.MakeTranslation(Cylinder_1, "l_x", "l_y", 0)
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
Regular_1D = Mesh_1.Segment()
Number_of_Segments_1 = Regular_1D.NumberOfSegments("l_x_p")
MEFISTO_2D = Mesh_1.Triangle(algo=smeshBuilder.MEFISTO)
status = Mesh_1.RemoveHypothesis(Number_of_Segments_1)
Number_of_Segments_2 = Regular_1D.NumberOfSegments("l_x_p")
status = Mesh_1.RemoveHypothesis(Regular_1D)
status = Mesh_1.RemoveHypothesis(Number_of_Segments_2)
NETGEN_3D = Mesh_1.Tetrahedron()
status = Mesh_1.RemoveHypothesis(MEFISTO_2D)
status = Mesh_1.RemoveHypothesis(NETGEN_3D)
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Simple_Parameters_1 = NETGEN_1D_2D_3D.Parameters(smeshBuilder.SIMPLE)
NETGEN_3D_Simple_Parameters_1.SetNumberOfSegments( "l_x_p" )
NETGEN_3D_Simple_Parameters_1.LengthFromEdges()
NETGEN_3D_Simple_Parameters_1.LengthFromFaces()
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0)
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/arani/input/assignment_cylinder_triang.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.FACE, 'Group_1_0' )
nbAdd = Group_1_0.Add( [ 13, 14, 15 ] )
[ Group_1_0 ] = Mesh_1.GetGroups()
Group_2_0 = Mesh_1.CreateEmptyGroup( SMESH.FACE, 'Group_1_1' )
nbAdd = Group_2_0.Add( [ 16, 17 ] )
[ Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Group_3_0 = Mesh_1.CreateEmptyGroup( SMESH.FACE, 'Group_1_2' )
nbAdd = Group_3_0.Add( [ 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40 ] )
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/arani/input/mesh_assignment_cylinder_triang.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/arani/input/mesh_assignment_cylinder_triang.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
Mesh_2 = smesh.Mesh(Translation_1)
Regular_1D_1 = Mesh_2.Segment()
Number_of_Segments_3 = Regular_1D_1.NumberOfSegments("l_x_p")
MEFISTO_2D_1 = Mesh_2.Triangle(algo=smeshBuilder.MEFISTO)
MG_Tetra = Mesh_2.Tetrahedron(algo=smeshBuilder.MG_Tetra)
status = Mesh_2.RemoveHypothesis(Number_of_Segments_3)
Number_of_Segments_4 = Regular_1D_1.NumberOfSegments(1)
status = Mesh_2.RemoveHypothesis(Number_of_Segments_4)
status = Mesh_2.RemoveHypothesis(MG_Tetra)
MG_Tetra_HPC = Mesh_2.Tetrahedron(algo=smeshBuilder.MG_Tetra_Parallel)
status = Mesh_2.RemoveHypothesis(MG_Tetra_HPC)
MG_Hybrid = Mesh_2.Tetrahedron(algo=smeshBuilder.HYBRID)
status = Mesh_2.RemoveHypothesis(MG_Hybrid)
NETGEN_3D_1 = Mesh_2.Tetrahedron()
Number_of_Segments_5 = Regular_1D_1.NumberOfSegments("l_x_p")
isDone = Mesh_2.Compute()
Mesh_2.ConvertToQuadratic(0)
Group_1_0_1 = Mesh_2.CreateEmptyGroup( SMESH.FACE, 'Group_1_0' )
nbAdd = Group_1_0_1.Add( [ 83, 84, 85, 86 ] )
[ Group_1_0_1 ] = Mesh_2.GetGroups()
Group_2_0_1 = Mesh_2.CreateEmptyGroup( SMESH.FACE, 'Group_2_0' )
nbAdd = Group_2_0_1.Add( [ 79, 80, 81, 82 ] )
[ Group_1_0_1, Group_2_0_1 ] = Mesh_2.GetGroups()
Group_3_0_1 = Mesh_2.CreateEmptyGroup( SMESH.FACE, 'Group_3_0' )
nbAdd = Group_3_0_1.Add( [ 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78 ] )
[ Group_1_0_1, Group_2_0_1, Group_3_0_1 ] = Mesh_2.GetGroups()
Group_2_0.SetName( 'Group_2_0' )
Group_3_0.SetName( 'Group_3_0' )
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/arani/input/mesh_assignment_cylinder_triang.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
smesh.SetName(Mesh_2, 'Mesh_2')
try:
  Mesh_2.ExportMED(r'/home/student/software/femus/applications/2021_fall/arani/input/mesh_assignment_cylinder_triang2.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(NETGEN_3D.GetAlgorithm(), 'NETGEN 3D')
smesh.SetName(MEFISTO_2D.GetAlgorithm(), 'MEFISTO_2D')
smesh.SetName(MG_Tetra.GetAlgorithm(), 'MG-Tetra')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(MG_Hybrid.GetAlgorithm(), 'MG-Hybrid')
smesh.SetName(MG_Tetra_HPC.GetAlgorithm(), 'MG-Tetra_HPC')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(NETGEN_3D_Simple_Parameters_1, 'NETGEN 3D Simple Parameters_1')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')
smesh.SetName(Number_of_Segments_5, 'Number of Segments_5')
smesh.SetName(Group_3_0_1, 'Group_3_0')
smesh.SetName(Group_2_0_1, 'Group_2_0')
smesh.SetName(Group_1_0_1, 'Group_1_0')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
