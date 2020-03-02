#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.3.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/applications/OptimalControl/pddm/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("lx", 1)
notebook.set("ly", 1)
notebook.set("lz", 1)
notebook.set("nx", 2)
notebook.set("ny", 3)
notebook.set("nz", 4)
####################################################
##        End of NoteBook variables section       ##
####################################################
###
### SHAPER component
###

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()
model.end()
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
Box_1 = geompy.MakeBoxDXDYDZ("lx", "ly", "lz")
[Face_1,Face_2,Face_3,Face_4,Face_5,Face_6] = geompy.ExtractShapes(Box_1, geompy.ShapeType["FACE"], True)
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Face_2, geompy.ShapeType["EDGE"], True)
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_3, geompy.ShapeType["EDGE"], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudyInFather( Box_1, Face_1, 'Face_1' )
geompy.addToStudyInFather( Box_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Box_1, Face_3, 'Face_3' )
geompy.addToStudyInFather( Box_1, Face_4, 'Face_4' )
geompy.addToStudyInFather( Box_1, Face_5, 'Face_5' )
geompy.addToStudyInFather( Box_1, Face_6, 'Face_6' )
geompy.addToStudyInFather( Face_3, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_3, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_3, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Face_3, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( Face_2, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_2, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_2, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Face_2, Edge_8, 'Edge_8' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Box_1)
Regular_1D = Mesh_1.Segment()
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Hexa_3D = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
Regular_1D_1 = Mesh_1.Segment(geom=Edge_2)
Regular_1D_2 = Mesh_1.Segment(geom=Edge_3)
#status = Mesh_1.AddHypothesis(smeshObj_1,Edge_3) ### smeshObj_1 has not been yet created
Regular_1D_3 = Mesh_1.Segment(geom=Edge_5)
#status = Mesh_1.AddHypothesis(smeshObj_2,Edge_3) ### smeshObj_2 has not been yet created
status = Mesh_1.RemoveHypothesis(Regular_1D)
status = Mesh_1.RemoveHypothesis(Quadrangle_2D)
status = Mesh_1.RemoveHypothesis(Hexa_3D)
Cartesian_3D = Mesh_1.BodyFitted()
#Mesh_1.GetMesh().RemoveSubMesh( smeshObj_3 ) ### smeshObj_3 has not been yet created
status = Mesh_1.RemoveHypothesis(Cartesian_3D)
Hexa_3D_1 = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
Regular_1D_4 = Mesh_1.Segment(geom=Edge_1)
Regular_1D_5 = Mesh_1.Segment()
Quadrangle_2D_1 = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
status = Mesh_1.RemoveHypothesis(Hexa_3D)
Hexa_3D_2 = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
status = Mesh_1.RemoveHypothesis(Hexa_3D)
Hexa_3D_3 = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
Number_of_Segments_7 = Regular_1D_1.NumberOfSegments("nx")
Propagation_of_1D_Hyp = Regular_1D_1.Propagation()
Number_of_Segments_8 = Regular_1D_3.NumberOfSegments("nz")
Propagation_of_1D_Hyp_1 = Regular_1D_3.Propagation()
Number_of_Segments_9 = Regular_1D_4.NumberOfSegments("ny")
Propagation_of_1D_Hyp_2 = Regular_1D_4.Propagation()
isDone = Mesh_1.Compute()
Sub_mesh_1 = Regular_1D_1.GetSubMesh()
Sub_mesh_3 = Regular_1D_3.GetSubMesh()
Sub_mesh_2 = Regular_1D_4.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_5')
smesh.SetName(Number_of_Segments_8, 'Number of Segments_8')
smesh.SetName(Propagation_of_1D_Hyp_1, 'Propagation of 1D Hyp. on Opposite Edges_6')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Number_of_Segments_9, 'Number of Segments_9')
smesh.SetName(Propagation_of_1D_Hyp_2, 'Propagation of 1D Hyp. on Opposite Edges_7')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Cartesian_3D.GetAlgorithm(), 'Cartesian_3D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Number_of_Segments_7, 'Number of Segments_7')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
