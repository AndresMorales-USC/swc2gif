# swc2gif
Converts SWC(s) (and voltage and other data) into a VTK model(s) and then into an animated GIF(s)

# I recommend setting up a virtual environment using conda (not entirely necissary)
conda create -n swc2gif_env 
# Use instead, if not using Guanbo's neuron install steps (which includes installing python 3.6)
conda create -n swc2gif_env python=3.6 anaconda

# Activate the virtual environment
source activate swc2gif_env

# INSTALLATION Packages
# Install vtk packages using conda
conda install -c anaconda vtk

# Install mayavi and prerequisites
conda install -c conda-forge pyside=1.2.4
pip install numpy scipy pygments imageio
pip install tqdm envisage PySide mayavi==4.5.0
