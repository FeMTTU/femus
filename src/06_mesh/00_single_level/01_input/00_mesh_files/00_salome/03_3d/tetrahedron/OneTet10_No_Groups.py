#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/src/06_mesh/00_single_level/01_input/00_mesh_files/00_salome/03_3d/tetrahedron')

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
Vertex_1 = geompy.MakeVertex(0, 1, 0)
Vertex_2 = geompy.MakeVertex(1, 0, 0)
Vertex_3 = geompy.MakeVertex(0, 0, 1)
Line_1 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_1)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_1)
Line_3 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_2)
Vertex_4 = geompy.MakeVertex(0, 0, 0)
Line_4 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_4)
Line_5 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_4)
Line_5_vertex_3 = geompy.GetSubShape(Line_5, [3])
Line_3_vertex_3 = geompy.GetSubShape(Line_3, [3])
Line_6 = geompy.MakeLineTwoPnt(Line_5_vertex_3, Line_3_vertex_3)
Face_1 = geompy.MakeFaceWires([Line_3, Line_4, Line_6], 1)
Face_2 = geompy.MakeFaceWires([Line_1, Line_2, Line_3], 1)
Face_3 = geompy.MakeFaceWires([Line_1, Line_4, Line_5], 1)
Face_4 = geompy.MakeFaceWires([Line_2, Line_5, Line_6], 1)
Shell_1 = geompy.MakeShell([Face_1, Face_2, Face_3, Face_4])
Solid_1 = geompy.MakeSolid([Shell_1])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudyInFather( Line_5, Line_5_vertex_3, 'Line_5:vertex_3' )
geompy.addToStudyInFather( Line_3, Line_3_vertex_3, 'Line_3:vertex_3' )
geompy.addToStudy( Line_6, 'Line_6' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Face_2, 'Face_2' )
geompy.addToStudy( Face_3, 'Face_3' )
geompy.addToStudy( Face_4, 'Face_4' )
geompy.addToStudy( Shell_1, 'Shell_1' )
geompy.addToStudy( Solid_1, 'Solid_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Solid_1,'Mesh_1')
NETGEN_2D3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Simple_Parameters_1 = NETGEN_2D3D.Parameters(smeshBuilder.SIMPLE)
NETGEN_3D_Simple_Parameters_1.SetNumberOfSegments( 1 )
NETGEN_3D_Simple_Parameters_1.LengthFromEdges()
NETGEN_3D_Simple_Parameters_1.LengthFromFaces()
isDone = Mesh_1.Compute()
Mesh_1.ConvertToQuadratic(1)
Mesh_1.ConvertToQuadratic(1, Mesh_1,True)
Mesh_1.ConvertToQuadratic(1)
smesh.SetName(Mesh_1, 'Mesh_1')


## Set names of Mesh objects
smesh.SetName(NETGEN_2D3D.GetAlgorithm(), 'NETGEN_2D3D')
smesh.SetName(NETGEN_3D_Simple_Parameters_1, 'NETGEN 3D Simple Parameters_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
