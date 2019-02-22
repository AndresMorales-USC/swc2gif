# swc2gif
Converts SWC(s) (and voltage and other data) into a VTK model(s) and then into an animated GIF(s)


# Install mayavi and other prerequisites
    (I recommend setting up a virtual environment using conda, though not entirely necissary)
conda create -n swc2gif_env python=3.6 anaconda

source activate swc2gif_env

pip install mayavi

pip install --upgrade tqdm numpy scipy imageio


# Functions
    function import code
    func_name(required_argument(s), [optional_argument], [optional_keyword_arguments=o_k_a])
    function description

from swc2gif.mapswc import mapswc

mapswc(swc_path_file(s), data_path_file, coordinate_path_file,
    [save_path=swc_path, vtk_name='model.vtk', datatitle='vdata',
    datadelimiter=' ', coordsdelimiter='	'])
    
Takes morphology file(s) and maps nearest neighboring data to each
    section. Data for each swc is saved into a separate associated
    .v file and the file path list is returned.

from swc2gif.swc2vtk import swc2vtk

swc2vtk(swc_path_file(s), [data_path_file(s)],
    [save_path=swc_path, vtk_name='model.vtk', datatitle='vdata',
    datadelimiter=' ', coordsdelimiter='	', maxframes=float('Inf'),
    spherediv=6, cyldiv=8])
    
Takes swc morphology file(s) and converts them into vtk model files.
    Also, returns number of vtk's saved and the min and max data values
    of all frames of all vtk's.

from swc2gif.vtk2gif import vtk2gif

vtk2gif(vtk_path_file,
    [save_path=vtk_path, datatitle='vdata', size=(1000,1000),
    showaxes=True, colorBounds=(), time_path_file='', numvtks=1,
    selectframes=[] or '', maxframes=float('Inf'), angleradians=False,
    startelevation=0, startazimuth=0, stepelevation=0, stepazimuth=0])
    
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
