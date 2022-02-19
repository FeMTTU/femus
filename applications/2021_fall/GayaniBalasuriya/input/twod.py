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

###
### SHAPER component
###

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()
Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()
Point_2 = model.addPoint(Part_1_doc, 0, 0, 0)
Point_3 = model.addPoint(Part_1_doc, 1, 1, 0)
Sketch_1 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))
SketchLine_1 = Sketch_1.addLine(1, 0, 0, 0)
SketchProjection_1 = Sketch_1.addProjection(model.selection("VERTEX", "Point_1"), False)
SketchPoint_1 = SketchProjection_1.createdFeature()
SketchConstraintCoincidence_1 = Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchPoint_1.result())
SketchLine_2 = Sketch_1.addLine(0, 0, 0, 1)
SketchLine_3 = Sketch_1.addLine(0, 1, 1, 1)
SketchLine_4 = Sketch_1.addLine(1, 1, 1, 0)
SketchConstraintCoincidence_2 = Sketch_1.setCoincident(SketchLine_4.endPoint(), SketchLine_1.startPoint())
SketchConstraintCoincidence_3 = Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchLine_2.startPoint())
SketchConstraintCoincidence_4 = Sketch_1.setCoincident(SketchLine_2.endPoint(), SketchLine_3.startPoint())
SketchConstraintCoincidence_5 = Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchLine_4.startPoint())
SketchConstraintHorizontal_1 = Sketch_1.setHorizontal(SketchLine_1.result())
SketchConstraintVertical_1 = Sketch_1.setVertical(SketchLine_2.result())
SketchConstraintHorizontal_2 = Sketch_1.setHorizontal(SketchLine_3.result())
SketchConstraintVertical_2 = Sketch_1.setVertical(SketchLine_4.result())
SketchProjection_2 = Sketch_1.addProjection(model.selection("VERTEX", "Point_2"), False)
SketchPoint_2 = SketchProjection_2.createdFeature()
SketchConstraintCoincidence_6 = Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchPoint_2.result())
model.do()
Export_1 = model.exportToXAO(Part_1_doc, '/tmp/shaper_hyp5lmo6.xao', model.selection(), 'XAO')
Export_2 = model.exportToXAO(Part_1_doc, '/tmp/shaper_zim4olsq.xao', model.selection(), 'XAO')
Export_3 = model.exportToXAO(Part_1_doc, '/tmp/shaper_8m2nc5fq.xao', model.selection(), 'XAO')
Export_4 = model.exportToXAO(Part_1_doc, '/tmp/shaper_kg3o6323.xao', model.selection(), 'XAO')
Face_1 = model.addFace(Part_1_doc, [model.selection("FACE", "Sketch_1/Face-SketchLine_4r-SketchLine_3r-SketchLine_2r-SketchLine_1r")])
Export_5 = model.exportToXAO(Part_1_doc, '/tmp/shaper_1umb9sl4.xao', model.selection(), 'XAO')
Export_6 = model.exportToXAO(Part_1_doc, '/tmp/shaper_k7_v2yrx.xao', model.selection(), 'XAO')
Export_7 = model.exportToXAO(Part_1_doc, '/tmp/shaper_q68w58xt.xao', model.selection(), 'XAO')
Export_8 = model.exportToXAO(Part_1_doc, '/tmp/shaper_cslox0vm.xao', model.selection(), 'XAO')
Export_9 = model.exportToXAO(Part_1_doc, '/tmp/shaper_ejab4pdf.xao', model.selection("FACE", "Face_1_1"), 'XAO')
Axis_4 = model.addAxis(Part_1_doc, 1, 0, 0)
Rotation_1 = model.addRotation(Part_1_doc, [model.selection("FACE", "Face_1_1")], model.selection("EDGE", "all-in-Axis_1"), 90)
Export_10 = model.exportToXAO(Part_1_doc, '/tmp/shaper_6t42384k.xao', model.selection("FACE", "Rotation_1_1"), 'XAO')
model.do()
model.end()
###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

(imported, Face_1_1, [], [], []) = geompy.ImportXAO("/tmp/shaper_ejab4pdf.xao")
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_1_1, geompy.ShapeType["EDGE"], True)
(imported, Rotation_1_1, [], [], []) = geompy.ImportXAO("/tmp/shaper_6t42384k.xao")
[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Rotation_1_1, geompy.ShapeType["EDGE"], True)
geompy.addToStudy( Face_1_1, 'Face_1_1' )
geompy.addToStudyInFather( Face_1_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_1_1, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_1_1, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Face_1_1, Edge_4, 'Edge_4' )
geompy.addToStudy( Rotation_1_1, 'Rotation_1_1' )
geompy.addToStudyInFather( Rotation_1_1, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Rotation_1_1, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Rotation_1_1, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Rotation_1_1, Edge_8, 'Edge_8' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Rotation_1_1)
Regular_1D = Mesh_1.Segment(geom=Edge_5)
Sub_mesh_1 = Regular_1D.GetSubMesh()
Number_of_Segments_1 = Regular_1D.NumberOfSegments(2)
Propagation_of_1D_Hyp = Regular_1D.Propagation()
Regular_1D_1 = Mesh_1.Segment(geom=Edge_6)
Sub_mesh_2 = Regular_1D_1.GetSubMesh()
Number_of_Segments_2 = Regular_1D_1.NumberOfSegments(2)
Propagation_of_1D_Hyp_1 = Regular_1D_1.Propagation()
Regular_1D_2 = Mesh_1.Segment()
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_1.Compute()
Group_1_0 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Group_1_0' )
nbAdd = Group_1_0.AddFrom( Mesh_1.GetMesh() )
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'/home/gbornia/software/femus/applications/OptimalControl/pddm/input/Mesh_2_xy.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=0)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Mesh_1.ConvertToQuadratic(0, Sub_mesh_1)
Mesh_1.ConvertToQuadratic(0, Sub_mesh_2)
isDone = Mesh_1.Compute()
[ Group_1_0 ] = Mesh_1.GetGroups()
Mesh_1.ConvertToQuadratic(0, Mesh_1,True)
[ Group_1_0 ] = Mesh_1.GetGroups()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_1_0, 'Group_1_0')
smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
smesh.SetName(Propagation_of_1D_Hyp, 'Propagation of 1D Hyp. on Opposite Edges_1')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Propagation_of_1D_Hyp_1, 'Propagation of 1D Hyp. on Opposite Edges_2')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
