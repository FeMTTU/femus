#!/bin/bash


echo -e "==============================
This script will 
- create a \$HOME/software/ folder and download Paraview there (~300MB), 
- extract the compressed file,
- install a wrapper script in your \$HOME/bin/ folder. \n
This way you will be able to launch Paraview from any bash terminal just by typing \"paraview\" 
============================== \n"
################################################
############# PARAVIEW #########################
################################################

export SW_DIR=software
cd 
mkdir -p $SW_DIR
cd $SW_DIR

# Download and extract paraview

export PV_BINNAME=paraview

export PV_FILEPREFIX=ParaView-5.6.0-MPI-Linux-64bit
wget -O  ${PV_FILEPREFIX}.tar.gz "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.6&type=binary&os=Linux&downloadFile=ParaView-5.6.0-MPI-Linux-64bit.tar.gz"  # beware of symbols interpreted by the shell!, and beware of the meaning of "-O" which is more than what we use it for

tar xzvf ${PV_FILEPREFIX}.tar.gz 

# create the script in the bin/ directory
cd 
mkdir -p bin/
cd bin
# Here, remember to put SLASH DOLLARS where needed to escape!!!
cat > $PV_BINNAME << EOF
#!/bin/bash
export PV_PATH=\$HOME/$SW_DIR/$PV_FILEPREFIX
\$PV_PATH/bin/$PV_BINNAME
EOF

# make the script executable
chmod a+x $PV_BINNAME

