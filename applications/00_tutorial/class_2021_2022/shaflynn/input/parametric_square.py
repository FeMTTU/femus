#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/max/software/femus/applications/2021_fall/shaflynn/input')

####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("l_y", 1)
notebook.set("l_x", 2)
notebook.set("half_l_y", "0.5*l_y")
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

Face_1 = geompy.MakeFaceHW("l_y", "l_x", 1)
Face_2 = geompy.MakeFaceHW("l_x", "l_y", 1)
Translation_1 = geompy.MakeTranslation(Face_2, "l_x", 0, 0)
Translation_2 = geompy.MakeTranslation(Face_2, "l_x", "half_l_y", 0)
Translation_3 = geompy.MakeTranslation(Face_2, "l_x", "half_l_y", 0)
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Face_2, 'Face_2' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudy( Translation_2, 'Translation_2' )
geompy.addToStudy( Translation_3, 'Translation_3' )
geompy.hideInStudy(Face_1)
geompy.hideInStudy(Translation_1)
geompy.hideInStudy(Translation_2)


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
