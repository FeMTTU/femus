import os

import medcoupling as MC
import MEDLoader as ML

# med file from the SAMPLES directory
#samples_dir = os.getenv("DATA_DIR")
#filename = os.path.join(samples_dir, "MedFiles", "fra.med")
filename = os.path.join("/home/gbornia/software/femus/unittests/test_mesh_read_write/input/" , "turek_FSI1.med") #I need to provide the absolute path

# Import the med file
UMesh = ML.MEDFileUMesh.New(filename)

# Get the mesh at level 0 (i.e. the volume mesh)
meshLvl0 = UMesh.getMeshAtLevel(0)

# Creates the faces around each cell by descending connectivity
# The face mesh is at level -1
meshLvlm1,d0,d1,d2,d3=meshLvl0.buildDescendingConnectivity()

# Creates the edges around each cell by descending connectivity
# The edge mesh is at level -2
#meshLvlm2,d0,d1,d2,d3=meshLvlm1.buildDescendingConnectivity()

# Write the med file
# We can write all cells or only the cells at one level
meshMEDFile=ML.MEDFileUMesh.New()
# volume mesh
meshMEDFile.setMeshAtLevel(0,meshLvl0)
# face mesh
meshMEDFile.setMeshAtLevel(-1,meshLvlm1)
# edges mesh
#meshMEDFile.setMeshAtLevel(-2,meshLvlm2)
meshMEDFile.write("fra_edges.med", 2) # 2 stands for write from scratch

