#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/student/software/femus/applications/2021_fall/hdongamm/input')

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
Disk_1 = geompy.MakeDiskR(1, 1)
geomObj_1 = geompy.MakeMarker(0, 0, 0, 1, 0, 0, 0, 1, 0)
sk = geompy.Sketcher2D()
sk.addPoint(-0.500000, -0.500000)
sk.addSegmentAbsolute(-0.500000, 0.500000)
sk.addSegmentAbsolute(0.500000, 0.500000)
sk.addSegmentAbsolute(0.500000, -0.500000)
sk.close()
Sketch_1 = sk.wire(Disk_1)
Wire_1 = geompy.MakeWire([Sketch_1], 1e-07)
Vertex_1 = geompy.MakeVertex(0, 0, 0)
Vertex_2 = geompy.MakeVertex(0.5, 0.5, 0)
Vertex_3 = geompy.MakeVertex(-0.5, 0.5, 0)
Vertex_4 = geompy.MakeVertex(0.707106, 0.707106, 0)
Vertex_5 = geompy.MakeVertex(-0.707106, 0.707106, 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_5)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_4)
Face_1 = geompy.MakeFaceWires([Sketch_1], 1)
Cut_1 = geompy.MakeCutList(Disk_1, [Face_1], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Disk_1, 'Disk_1' )
geompy.addToStudy( Sketch_1, 'Sketch_1' )
geompy.addToStudy( Wire_1, 'Wire_1' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Cut_1, 'Cut_1' )


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
