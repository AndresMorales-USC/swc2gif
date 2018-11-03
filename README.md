# swc2gif
Convert SWC(s) (and voltage and other data) into a VTK model(s) and then into an animated GIF(s)

# I recommend setting up a virtual environment using conda (not entirely necissary)
conda create -n swc2gif_env python=2.7 anaconda
source activate swc2gif_env

# INSTALLATION prerequisites
# Install vtk packages using conda
conda install -c anaconda vtk

# Install tqdm and mayavi (mayavi requires envisage and PySide)
conda install -c conda-forge pyside=1.2.4 #linux-py3.6 install
pip install numpy scipy pygments imageio #linux-py3.6 install
pip install tqdm envisage PySide mayavi==4.5.0
