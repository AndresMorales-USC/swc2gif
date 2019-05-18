# swc2gif
Converts SWC(s) (and voltage and other data) into a VTK model(s) and then into an animated GIF(s)


# Install mayavi and other prerequisites
(I recommend setting up a virtual environment using conda, though not entirely necissary)
    
    conda create -n swc2gif_env python=3.6
    conda activate swc2gif_env
    pip install numpy vtk tqdm scipy imageio
    pip install mayavi
    pip install PyQt5

# Possible installation issues:
Conda environments fail to activate (CommandNotFoundError: No command 'conda conda'.) and possibly break conda command
- Solution: revert conda to older version (4.6.7 or 4.5.12) before creating environment

        conda install -n base conda==4.6.7

PyQt5 causing windows error about Bluetooth API
- Solution: Install older version of PyQt5 (5.9.2)

        pip install PyQt5==5.9.2


# Functions
    function import code
    func_name(required_argument(s), [optional_arguments], [optional_keyword_arguments=o_k_a])
    function description

from swc2gif.swc2vtk import swc2vtk
    
swc2vtk(swc_path_file(s), [data_path_file, coord_path_file],
    [
    time_path_file='', save_path=swc_path, vtk_name='model.vtk',
    datatitle='vdata', datadelimiter=' ', coordsdelimiter='	',
    fpvtk=float('Inf'), spherediv=6, cyldiv=8,
    invertSWC=((False, False, False),), scaleSWC=(1.0,),
    shiftSWC=((0.0,0.0,0.0),)
    ])

Takes swc morphology file(s) and converts them into vtk model files.
    Also, returns the min and max data values of all vtk's.

from swc2gif.vtk2gif import vtk2gif
    
vtk2gif(vtk_path_file,
    [
    save_path=vtk_path, datatitle='vdata', size=(1000,1000),
    showaxes=True, colorBounds=(), time_list=[], selectframes=[],
    fpvtk=0, totalframes=0, angleradians=False,
    startelevation=0, startazimuth=0, stepelevation=0, stepazimuth=0,
    bgcolor=(0.0,0.0,0.0)
    ])

Takes VTK file(s) containing vdata, where each vdata## is a different
    frame of the produced GIF(s).


# Example Code
(requires editing variables with *_path* to point to proper input files/folders)

    python -m swc2gif.example.py

or

    ipython -m swc2gif.example.py
or

    ipython swc2gif/example.py
or

    mayavi2 -x 'swc2gif/example.py'
