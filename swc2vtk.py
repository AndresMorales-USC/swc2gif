# -*- coding: utf-8 -*-
"""
Created on 2018/05/22

@author: Andres Morales

swc2vtk(swc_path_file(s), [data_path_file(s)],
    [
    save_path=swc_path, vtk_name='model.vtk', datatitle='vdata',
    datadelimiter=' ', coordsdelimiter='	', maxframes=float('Inf'),
    spherediv=6, cyldiv=8, invertSWC=[False, False, False],
    scaleSWC=1.0, shiftSWC=[0.0,0.0,0.0]
    ])

Takes swc morphology file(s) and converts them into vtk model files.
    Also, returns number of vtk's saved and the min and max data values
    of all frames of all vtk's.


Required Argument:
swc_path_file(s): A string containing the file path and file name,
    or tuple/list of strings containing the file paths and file names,
    of the ".swc" file(s).
    (Example 1: 'C:\\Retina simulation\\morphology\\606.swc'
     Example 2: ('C:\\Retina simulation\\morphology\\606.swc',
                'C:\\Retina simulation\\morphology\\1010.swc')
            or: ['C:\\Retina simulation\\morphology\\606.swc',
                'C:\\Retina simulation\\morphology\\1010.swc'])


Optional Arguments:
data_path_file(s): A string containing the file path and file name,
    or tuple/list of strings containing the file paths and file names,
    of the ".txt" data file(s). Each data file is information to be
    mapped to a specific swc (swc & data file lists must be of equal
    size and ordering).
    (Example 1: 'C:\\Retina simulation\\morphology\\voltage_data.txt'
     Example 2: ['C:\\Retina simulation\\morphology\\606_data.txt',
                 'C:\\Retina simulation\\morphology\\1010_data.txt'])


Optional Kew Word Arguments:
save_path: A string containing the file path of where to save generated
    files (alinged data files or ".vtk" files). Save_path defaults to 
    the same folder as where the ".swc" files are read.
    (Example: 'C:\\Retina simulation\\outputs')

vtk_name: A string containing the file name of the ".vtk" file to be
    generated. vtk_name, if not given, is set to 'model.vtk'.

datatitle: A string containing the name to label data within the ".vtk"
    file. It defaults to 'vdata'.

datadelimiter: A string used to split input from the data file. The
    default data delimiter is ' '.

maxframes: An integer for the maximum number of frames to
    be written to a single ".vtk". If the number of frames in the data
    file exceedes maxframes, multiple ".vtk"s will be written. Each vtk
    containing at most the frame quota and sequentially labeled with
    the prefix '_' and a number (starting at 0). The default value is
    actually a float for infinity.
    
spherediv: An integer for the number of radial sides to modeled spheres.
    The default value is 6.
    
cyldiv:  An integer for the number of radial sides to modeled cylinders.
    The default value is 8.

invertSWC: A list of 3 booleans for inverting the axis(es) of the swc
    coordinates. This is done before scaling and shifting of swc
    coordinates. Default is [False, False, False].

scaleSWC: A float for scaling the swc coordinates to fit the
    dimensions of the data files. It is done after inversion and before
    shifting of swc coordinates. Default is 1.0.
    
shiftSWC: A list of 3 floats for shifting the origin of the swc
    coordinates. This is done after scaling and inversion of swc
    coordinates. Default is [0.0, 0.0, 0.0].


Prerequisite Packages: tqdm
pip install tqdm

"""

import swc2gif.vtkgen as vtkgen
import os
from tqdm import tqdm

def swc2vtk(*args, **kwargs):
    # Check if swc and data file paths were passed through arg
    # If only 1, assume data was intentionally not included
    # If more than 2 will ignore extra args
    if len(args) == 0:
        print('ERROR: Not enough arguments [minimum: list of swc''s]')
    else: # SWC arguments are contained in args
        # Place file path strings into a list from whatever format it was submitted
        if isinstance(args[0], str):
            swc_list = [args[0]]
        elif isinstance(args[0], tuple):
            swc_list = list(args[0])
        elif isinstance(args[0], list):
            swc_list = args[0]
        else:
            print('Warning: Unknown input for swc file(s)')
        
        # Check if all SWC's have swc as file types
        if not checktype(swc_list, 'swc'):
            print('Warning: One or more SWC\'s do NOT have a \'.swc\' file type')
        
        data_file_list = []
        
        if len(args) > 1: # Data arguments are contained in args
            if isinstance(args[1], str):
                data_file_list = [args[1]]
            elif isinstance(args[1], tuple):
                data_file_list = list(args[1])
            elif isinstance(args[1], list):
                data_file_list = args[1]
            else:
                print('Warning: Unknown input for data file(s)')
            # Check if all data V's have v as file types
            if not checktype(data_file_list, 'v'):
                print('Warning: One or more data V\'s do NOT have a \'.v\' file type')
            # If data V's are used, check if number is appropriate for SWC's
            if len(data_file_list)>1 and len(data_file_list) != len(swc_list):
                print('Warning: Mismatch in number of SWC and data V files')
            
            if len(args) > 2: # Extra arguments are contained in args
                print('Warning: Extra arguments recieved beyond SWC and Data files')
        else: # Less than 2 arguments
            print('No data files. Converting swc''s only.')
                
        
    
    # Check if keys are defined in kwargs
    # Else set to defaults
    
    # Path to save VTK file default: '' (results in saving in same folder as swc's)
    save_path = kwargs.get('save_path', '')
    
    # VTK file name default: 'model.vtk'
    vtk_name = kwargs.get('vtk_name', 'model.vtk')
    
    # Title of data default: 'vdata'
    datatitle = kwargs.get('datatitle', 'vdata')
    
    # Data file default delimiter: ' '
    datadelimiter = kwargs.get('datadelimiter', ' ')
    
    # Coordinate file default delimiter: '	'
    coordsdelimiter = kwargs.get('coordsdelimiter', '	')
    
    # Maximum frames per vtk default: infinite
    maxframes = kwargs.get('maxframes', float('Inf'))

    # Default # of divisions of modeled spheres: 6
    spherediv = kwargs.get('spherediv', 6)
    
    # Default # of divisions of modeled cylinders: 8
    cyldiv = kwargs.get('cyldiv', 8)

    # invertSWC = [x, y, z] [False, False, False]:
    # Done before scale and shift of swc coordinates
    invertSWC = kwargs.get('invertSWC', [False, False, False])
    
    # scaleSWC = 1.0:
    # scale size (positions and radius) of swc
    # Done after invert and before shift of swc coordinates
    scaleSWC = kwargs.get('scaleSWC', 1.0)
    
    # shiftSWC = [x,y,z] [0.0, 0.0, 0.0]:
    # Done after invert and scale of swc coordinates
    shiftSWC = kwargs.get('shiftSWC', [0.0, 0.0, 0.0])
    
    # Initialize VTK and add SWC file(s)
    head, tail = os.path.split(swc_list[0])
    # If no save path was specified use the swc's path
    if len(save_path) == 0:
        vtk_path_file = os.path.join(head, vtk_name)
    else:
        vtk_path_file = os.path.join(save_path, vtk_name)
    vtk = vtkgen.VtkGenerator()
    for swc_path_file in swc_list:
        vtk.add_swc(swc_path_file, shift_x=shiftSWC[0], shift_y=shiftSWC[1], shift_z=shiftSWC[2],
                                    inv_x=invertSWC[0], inv_y=invertSWC[1], inv_z=invertSWC[2],
                                    scale_factor=scaleSWC)
    
    
    # Check if data needs to be added to VTK
    if len(data_file_list) > 0:
        # Add all data to VTK and write it to the same folder as the first SWC
        for data_path_file in data_file_list:
            vtk.add_timeddatafile(data_path_file)
        numvtks, minV, maxV = vtk.write_vtk(vtk_path_file, datatitle=datatitle, delimiter=datadelimiter, maxframes=maxframes, sphere_div=spherediv, cyl_div=cyldiv)
        
    else: # No data to add, just write the VTK to the same folder as the first SWC
        numvtks, minV, maxV = vtk.write_vtk(vtk_path_file, sphere_div=sphere_div, cyl_div=cyl_div)
    
    # Save split data
    print('\nNumber of VTK''s: '+str(numvtks))
    print('Data Bounds [Min, Max]: ['+str(minV)+', '+str(maxV)+']')
    # If no save path was specified use the swc's path
    if len(save_path) == 0:
        vtksplitdata_path_file = os.path.join(head, vtk_name[:-4]+'_splitdata.txt')
    else:
        vtksplitdata_path_file = os.path.join(save_path, vtk_name[:-4]+'_splitdata.txt')
    with open(vtksplitdata_path_file, 'w') as f:
        f.write('Number of VTK''s: '+str(numvtks)+'\n')
        f.write('Data Bounds [Min, Max]: ['+str(minV)+', '+str(maxV)+']')
    return [numvtks, minV, maxV]





def checktype(file_list, file_type):
    # Returns true if all files in file_list are of type file_type
    file_type_check = True
    chars = len(file_type)
    for path_file in file_list:
        if path_file[-chars:] != file_type:
            file_type_check = False
            break
    return file_type_check
