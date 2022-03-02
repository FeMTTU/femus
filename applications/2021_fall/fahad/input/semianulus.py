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
notebook.set("n_x", 6)
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

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Cut_2)
QuadFromMedialAxis_1D2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUAD_MA_PROJ)
Number_of_Layers_1 = QuadFromMedialAxis_1D2D.NumberOfLayers("n_x")
isDone = Mesh_1.Compute()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/fine_mesh_semiannulus.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/student/software/femus/applications/2021_fall/fahad/input/Mesh_fine_semiannulus.med',auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(QuadFromMedialAxis_1D2D.GetAlgorithm(), 'QuadFromMedialAxis_1D2D')
smesh.SetName(Number_of_Layers_1, 'Number of Layers_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
