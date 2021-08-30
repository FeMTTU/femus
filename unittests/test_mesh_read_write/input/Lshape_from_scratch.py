#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.6.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/gbornia/software/femus/unittests/test_mesh_read_write/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("l_x", 1)
notebook.set("l_y", 1)
notebook.set("leg_distance_from_x_axis", 0.5)
notebook.set("leg_distance_from_y_axis", 0.5)
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

geomObj_1 = geompy.MakeVertex(0, 0, 0)
geomObj_2 = geompy.MakeVectorDXDYDZ(1, 0, 0)
geomObj_3 = geompy.MakeVectorDXDYDZ(0, 1, 0)
geomObj_4 = geompy.MakeVectorDXDYDZ(0, 0, 1)
geomObj_5 = geompy.MakeVertex(0, 0, 0)
geomObj_6 = geompy.MakeVertex(1, 0, 0)
geomObj_7 = geompy.MakeVertex(1, 0.5, 0)
geomObj_8 = geompy.MakeVertex(0.5, 1, 0)
geomObj_9 = geompy.MakeVertex(0.5, 0.5, 0)
geomObj_10 = geompy.MakeVertex(0, 1, 0)
geomObj_11 = geompy.MakeVertex(0.5, 0, 0)
geomObj_12 = geompy.MakeVertex(0, 0.5, 0)
geomObj_13 = geompy.MakeLineTwoPnt(geomObj_5, geomObj_11)
geomObj_14 = geompy.MakeLineTwoPnt(geomObj_11, geomObj_9)
geomObj_15 = geompy.MakeLineTwoPnt(geomObj_9, geomObj_12)
geomObj_16 = geompy.GetSubShape(geomObj_13, [2])
geomObj_17 = geompy.MakeLineTwoPnt(geomObj_12, geomObj_16)
geomObj_18 = geompy.MakeLineTwoPnt(geomObj_11, geomObj_6)
geomObj_19 = geompy.MakeLineTwoPnt(geomObj_6, geomObj_7)
geomObj_20 = geompy.MakeLineTwoPnt(geomObj_7, geomObj_9)
geomObj_21 = geompy.GetSubShape(geomObj_18, [2])
geomObj_22 = geompy.MakeLineTwoPnt(geomObj_9, geomObj_21)
geomObj_23 = geompy.MakeLineTwoPnt(geomObj_12, geomObj_9)
geomObj_24 = geompy.MakeLineTwoPnt(geomObj_9, geomObj_8)
geomObj_25 = geompy.GetSubShape(geomObj_24, [3])
geomObj_26 = geompy.MakeLineTwoPnt(geomObj_25, geomObj_10)
geomObj_27 = geompy.MakeLineTwoPnt(geomObj_10, geomObj_12)
geomObj_28 = geompy.MakeFaceWires([geomObj_13, geomObj_14, geomObj_15, geomObj_17], 1)
geomObj_29 = geompy.MakeFaceWires([geomObj_18, geomObj_19, geomObj_20, geomObj_22], 1)
geomObj_30 = geompy.MakeFaceWires([geomObj_23, geomObj_24, geomObj_26, geomObj_27], 1)
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Vertex_1 = geompy.MakeVertex(0, 0, 0)
Vertex_2 = geompy.MakeVertex("leg_distance_from_y_axis", 0, 0)
Vertex_3 = geompy.MakeVertex("leg_distance_from_y_axis", "leg_distance_from_x_axis", 0)
Vertex_4 = geompy.MakeVertex(0, "leg_distance_from_x_axis", 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_3 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_4)
Line_4 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_1)
Face_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_3, Line_4], 1)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudy( Face_1, 'Face_1' )


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
