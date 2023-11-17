

salome start -t &

for file in `find . -name *.py`; do echo $file "==================="; salome shell $file; done
