 
for file in `find . -iname *.py`; do echo $file "==================="; salome_9.11.0 shell $file; done
