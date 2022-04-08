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
notebook.set("n_seg", 6)
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
Regular_1D = Mesh_1.Segment()
Number_of_Segments_1 = Regular_1D.NumberOfSegments("n_seg")
NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_2D)
NETGEN_3D = Mesh_1.Tetrahedron()
isDone = Mesh_1.Compute()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(NETGEN_3D.GetAlgorithm(), 'NETGEN 3D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
