#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/student/software/femus/applications/2021_fall/fahad/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("r_in", 0.5)
notebook.set("r_out", 1)
notebook.set("n_r", 1)
notebook.set("n_theta", 6)
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
Disk_1 = geompy.MakeDiskR("r_in", 1)
Disk_2 = geompy.MakeDiskR("r_out", 1)
Cut_1 = geompy.MakeCutList(Disk_2, [Disk_1], True)
Face_1 = geompy.MakeFaceHW(2, 2, 1)
Translation_1 = geompy.MakeTranslation(Face_1, -1, 0, 0)
Cut_2 = geompy.MakeCutList(Cut_1, [Translation_1], True)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Cut_2, geompy.ShapeType["EDGE"], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Disk_1, 'Disk_1' )
geompy.addToStudy( Disk_2, 'Disk_2' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudy( Cut_2, 'Cut_2' )
geompy.addToStudyInFather( Cut_2, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Cut_2, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Cut_2, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Cut_2, Edge_4, 'Edge_4' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Cut_2)
QuadFromMedialAxis_1D2D = smesh.CreateHypothesis('QuadFromMedialAxis_1D2D')
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/fahad/input/semiannulus_mar2.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Distribution_of_Layers_1 = smesh.CreateHypothesis('LayerDistribution2D')
GeometricProgression_Distribution = smesh.CreateHypothesis('GeometricProgression')
Distribution_of_Layers_1.SetLayerDistribution( GeometricProgression_Distribution )
Number_of_Layers_1 = smesh.CreateHypothesis('NumberOfLayers2D')
Number_of_Layers_1.SetNumberOfLayers( "n_r" )
Regular_1D = Mesh_1.Segment(geom=Edge_2)
Number_of_Segments_1 = Regular_1D.NumberOfSegments("n_r",None,[])
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Regular_1D_1 = Mesh_1.Segment(geom=Edge_3)
Number_of_Segments_2 = Regular_1D_1.NumberOfSegments("n_theta",None,[])
status = Mesh_1.AddHypothesis(Propagation_of_1D_Hyp,Edge_3)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Regular_1D_2 = Mesh_1.Segment()
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
Group_1_0 = Mesh_1.GroupOnGeom(Edge_2,'Group_1_0',SMESH.EDGE)
[ Group_1_0 ] = Mesh_1.GetGroups()
Group_2_0 = Mesh_1.GroupOnGeom(Edge_1,'Group_2_0',SMESH.EDGE)
[ Group_1_0, Group_2_0 ] = Mesh_1.GetGroups()
Group_3_0 = Mesh_1.GroupOnGeom(Edge_3,'Group_3_0',SMESH.EDGE)
[ Group_1_0, Group_2_0, Group_3_0 ] = Mesh_1.GetGroups()
Group_4_0 = Mesh_1.GroupOnGeom(Edge_4,'Group_4_0',SMESH.EDGE)
[ Group_1_0, Group_2_0, Group_3_0, Group_4_0 ] = Mesh_1.GetGroups()
Sub_mesh_2 = Mesh_1.GetSubMesh( Edge_2, 'Sub-mesh_2' )
Sub_mesh_3 = Mesh_1.GetSubMesh( Edge_3, 'Sub-mesh_3' )


## Set names of Mesh objects
smesh.SetName(QuadFromMedialAxis_1D2D, 'QuadFromMedialAxis_1D2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Group_2_0, 'Group_2_0')
smesh.SetName(Group_3_0, 'Group_3_0')
smesh.SetName(Group_4_0, 'Group_4_0')
smesh.SetName(Number_of_Layers_1, 'Number of Layers_1')
smesh.SetName(Distribution_of_Layers_1, 'Distribution of Layers_1')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
