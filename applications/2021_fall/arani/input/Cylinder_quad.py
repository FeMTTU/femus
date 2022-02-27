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
notebook.set("l_x_p", 1)
notebook.set("l_y_p", 1)
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
O_1 = geompy.MakeVertex(0, 0, 0)
OX_1 = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY_1 = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ_1 = geompy.MakeVectorDXDYDZ(0, 0, 1)
Divided_Disk_1 = geompy.MakeDividedDisk(1, 1, GEOM.SQUARE)
Rotation_1 = geompy.MakeRotation(Divided_Disk_1, OZ_1, 45*math.pi/180.0)
Translation_1 = geompy.MakeTranslation(Rotation_1, "l_x", "l_y", 0)
[Face_1,Face_2,Face_3,Face_4,Face_5] = geompy.ExtractShapes(Translation_1, geompy.ShapeType["FACE"], True)

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Number_of_Segments_1 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_1.SetNumberOfSegments( 4 )
Mesh_2 = smesh.Mesh(Face_1)
Regular_1D = Mesh_2.Segment()
Quadrangle_2D = Mesh_2.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_2.Compute()
Mesh_2.ConvertToQuadratic(0, Mesh_2,True)
Number_of_Segments_2 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_2.SetNumberOfSegments( 4 )
Mesh_1 = smesh.Mesh(Face_2)
Regular_1D_1 = Mesh_1.Segment()
Quadrangle_2D_1 = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_1.Compute()
Number_of_Segments_3 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_3.SetNumberOfSegments( "l_x" )
Mesh_3 = smesh.Mesh(Face_3)
Regular_1D_2 = Mesh_3.Segment()
Quadrangle_2D_2 = Mesh_3.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_3.Compute()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Mesh_3.ConvertToQuadratic(0, Mesh_3,True)
Number_of_Segments_4 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_4.SetNumberOfSegments( 4 )
isDone = Mesh_3.Compute()
Mesh_3.ConvertToQuadratic(0, Mesh_3,True)
Number_of_Segments_5 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_5.SetNumberOfSegments( 4 )
Mesh_4 = smesh.Mesh(Face_4)
Regular_1D_3 = Mesh_4.Segment()
Quadrangle_2D_3 = Mesh_4.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_4.Compute()
Mesh_4.ConvertToQuadratic(0, Mesh_4,True)
Number_of_Segments_6 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_6.SetNumberOfSegments( 4 )
Mesh_5 = smesh.Mesh(Face_5)
Regular_1D_4 = Mesh_5.Segment()
Quadrangle_2D_4 = Mesh_5.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_5.Compute()
Mesh_5.ConvertToQuadratic(0, Mesh_5,True)
Number_of_Segments_7 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_7.SetNumberOfSegments( "l_x" )
Number_of_Segments_8 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_8.SetNumberOfSegments( "l_x" )
Number_of_Segments_9 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_9.SetNumberOfSegments( "l_x" )
Number_of_Segments_10 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_10.SetNumberOfSegments( "l_x" )
Number_of_Segments_11 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_11.SetNumberOfSegments( "l_x" )
Number_of_Segments_12 = Regular_1D.NumberOfSegments("l_x_p")
isDone = Mesh_2.Compute()
Number_of_Segments_13 = Regular_1D_1.NumberOfSegments("l_x_p")
Number_of_Segments_14 = Regular_1D_2.NumberOfSegments("l_x_p")
Number_of_Segments_15 = Regular_1D_3.NumberOfSegments("l_x_p")
Number_of_Segments_16 = Regular_1D_4.NumberOfSegments("l_x_p")
isDone = Mesh_1.Compute()
isDone = Mesh_3.Compute()
isDone = Mesh_4.Compute()
isDone = Mesh_5.Compute()
Mesh_2.ConvertToQuadratic(0, Mesh_2,True)
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Mesh_3.ConvertToQuadratic(0, Mesh_3,True)
Mesh_4.ConvertToQuadratic(0, Mesh_4,True)
Mesh_5.ConvertToQuadratic(0, Mesh_5,True)
Mesh_6 = smesh.Concatenate( [ Mesh_2.GetMesh(), Mesh_1.GetMesh(), Mesh_3.GetMesh(), Mesh_4.GetMesh(), Mesh_5.GetMesh() ], 1, 1, 1e-05, False )
Mesh_6.ExtrusionSweepObjects( [ Mesh_6 ], [ Mesh_6 ], [ Mesh_6 ], [ 0, 0, 2 ], 1, 1, [  ], 0, [  ], [  ], 0 )
smesh.SetName(Mesh_6, 'Compound_Mesh_1')
smesh.SetName(Mesh_6, 'Mesh_6')
Group_1_0 = Mesh_6.CreateEmptyGroup( SMESH.FACE, 'Group_1_0' )
nbAdd = Group_1_0.Add( [ 5, 10, 15, 20, 25 ] )
[ Group_1_0 ] = Mesh_6.GetGroups()
Group_2_0 = Mesh_6.CreateEmptyGroup( SMESH.FACE, 'Group_2_0' )
nbAdd = Group_2_0.Add( [ 72, 73, 74, 75, 76 ] )
[ Group_1_0, Group_2_0 ] = Mesh_6.GetGroups()
Group_3_0 = Mesh_6.CreateEmptyGroup( SMESH.FACE, 'Group_3_0' )
nbAdd = Group_3_0.Add( [ 12, 21, 31, 33 ] )
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_6.GetGroups()
smesh.SetName(Mesh_6, 'Mesh_6')
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_6.GetGroups()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Number_of_Segments_3, 'Number of Segments_3')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Segments_6, 'Number of Segments_6')
smesh.SetName(Number_of_Segments_7, 'Number of Segments_7')
smesh.SetName(Number_of_Segments_4, 'Number of Segments_4')
smesh.SetName(Number_of_Segments_5, 'Number of Segments_5')
smesh.SetName(Number_of_Segments_8, 'Number of Segments_8')
smesh.SetName(Number_of_Segments_9, 'Number of Segments_9')
smesh.SetName(Mesh_6.GetMesh(), 'Mesh_6')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')
smesh.SetName(Mesh_3.GetMesh(), 'Mesh_3')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Mesh_5.GetMesh(), 'Mesh_5')
smesh.SetName(Mesh_4.GetMesh(), 'Mesh_4')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Number_of_Segments_13, 'Number of Segments_13')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Number_of_Segments_12, 'Number of Segments_12')
smesh.SetName(Number_of_Segments_11, 'Number of Segments_11')
smesh.SetName(Number_of_Segments_10, 'Number of Segments_10')
smesh.SetName(Number_of_Segments_16, 'Number of Segments_16')
smesh.SetName(Number_of_Segments_15, 'Number of Segments_15')
smesh.SetName(Number_of_Segments_14, 'Number of Segments_14')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
