# -*- coding: utf-8 -*-
"""
Created on 2019/01/05

@author: Andres Morales

mapswc(swc_path_file(s), data_path_file, coordinate_path_file,
    [save_path=swc_path, vtk_name='model.vtk', datatitle='vdata',
    datadelimiter=' ', coordsdelimiter='	'])

Takes morphology file(s) and maps nearest neighboring data to each
    section. Data for each swc is saved into a separate associated
    .v file and the file path list is returned.


Required Argument:
swc_path_file(s): A string containing the file path and file name,
    or tuple/list of strings containing the file paths and file names,
    of the ".swc" file(s).
    (Example: 'C:\\Retina simulation\\morphology\\606.swc'
          or: ('C:\\Retina simulation\\morphology\\606.swc',
               'C:\\Retina simulation\\morphology\\1010.swc')
          or: ['C:\\Retina simulation\\morphology\\606.swc',
               'C:\\Retina simulation\\morphology\\1010.swc'])

data_path_file: A string containing the file path and file name of the
    ".v" data file(s). Data must be in a single file where each line is 
    a different time step.
    (Example: 'C:\\Retina simulation\\morphology\\voltage_data.v')

coordinate_path_file(s): A string containing the file path and file
    name of the ".txt" coordinate file. Each line of the file is 3
    numbers separated by some delimiter giving a single location
    cooridnate
    (Example: 'C:\\Retina simulation\\morphology\\coords_data.txt')


Optional Kew Word Arguments:
save_path: A string containing the file path of where to save generated
    files. Save_path defaults to the same folder as where the ".swc"
    files are read.
    (Example: 'C:\\Retina simulation\\outputs')

datadelimiter: A string used to split input from the data file. The
    default data delimiter is ' '.

coordsdelimiter: A string used to split input from the coordinate file.
    The default coordinate data delimiter is '	'.


Prerequisite Packages: tqdm and scypi
pip install tqdm scypi

"""

import swc2gif.vtkgen as vtkgen
from swc2gif.vtkgen.swc import Swc
import os
from tqdm import tqdm
from scipy import spatial

def mapswc(*args, **kwargs):
    #TODO: Check if correct
    # Check if swc, data, and coordinate file paths were passed through arg
    # If more than 3 will ignore extra args
    if len(args) < 3:
        print('ERROR: Not enough arguments [minimum: morphology.swc, data.v, coordinates.txt')
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
        # Check if all SWC's have .swc as file types
        if not checktype(swc_list, 'swc'):
            print('Warning: One or more SWC\'s do NOT have a \'.swc\' file type')
        
        
        # Get data file path
        if isinstance(args[1], str):
            data_path_file = args[1]
        else:
            print('Warning: Unknown input for data file(s)')
        # Check if data V's have .v as file types
        if not checktype([data_path_file], 'v'):
            print('Warning: One or more data V\'s do NOT have a \'.v\' file type')
        
        # Get coordinate file path
        if isinstance(args[2], str):
            coords_path_file = args[2]
        else:
            print('Warning: Unknown input for coordinate file')
        # Check if all coordinate TXT's have txt as file types
        if not checktype([coords_path_file], 'txt'):
            print('Warning: One or more coordinate TXT\'s do NOT have a \'.txt\' file type')
            
        
        if len(args) > 3: # Extra arguments are contained in args
            print('Warning: Extra arguments recieved beyond SWC, Data, and Coordinate files')
    
    # Check if keys are defined in kwargs
    # Else set to defaults
    
    # Path to save VTK file default: '' (results in saving in same folder as swc's)
    save_path = kwargs.get('save_path', '')
    
    # Data file default delimiter: ' '
    datadelimiter = kwargs.get('datadelimiter', ' ')
    
    # Coordinate file default delimiter: '	'
    coordsdelimiter = kwargs.get('coordsdelimiter', '	')
    

    # Convert all data using SWC's and coordinate TXT's

    # Extract data coordinates from coords file
    tempHead, coords_file = os.path.split(coords_path_file)
    data_coords = []
    with open(coords_path_file, 'r') as f:
        # Read each line of the file as a string, use delimiter to split
        #     into a list of strings, convert to a list of numbers, and 
        #     finally append single postion list to list of all positions
        for line in tqdm(f, desc='Reading coordinate file'):
            str_list = line.rstrip().split(coordsdelimiter)
            float_list = []
            for val in str_list:
                try:
                    float_list.append(float(val))
                except:
                    print(val)
            data_coords.append(float_list)

    # Generate KDTree
    print('Generating KD Tree')
    tree = spatial.cKDTree(data_coords)

    aligned_data_file_list = []
    timesteps=0
    # Loop through SWC files for conversion
    for file_idx, swc_path_file in enumerate(tqdm(swc_list, desc='Aligning Data to SWC\'s')):
        tempHead, swc_file = os.path.split(swc_path_file)
        swc = Swc(swc_path_file)
        
        # Iterate through each swc compartment and generate reference list
        compartment_coords_idx = []
        for swc_compartment in swc.data.values():
            # Find index of nearest neighbor in coordinates and append to list
            distance, idx = tree.query(swc_compartment['pos'])
            compartment_coords_idx.append(idx)
        
        # Generate aligned data
        # If no save path was specified use the swc's path
        if len(save_path) == 0:
            aligned_data_path_file = swc_path_file[:-4]+'_aligned_data.v'
        else:
            aligned_data_path_file = os.path.join(save_path, swc_file[:-4]+'_aligned_data.v')
        with open(aligned_data_path_file, 'w') as adf:
            # Iterate through data file where each line is a different time step
            with open(data_path_file, 'r') as df:
                # Loop through all time steps
                linecount = 0
                for line in tqdm(df, total=timesteps, desc='Generating Aligned Data for '+swc_file):
                    aligned_list = []
                    #aligned_line = ''
                    original_str = line.rstrip().split(datadelimiter)
                    # Loop through each compartment and add nearest neighbor data value
                    for idx in compartment_coords_idx:
                        aligned_list.append(original_str[idx])
                        #aligned_line += original_str[idx]+datadelimiter
                    #aligned_line += '\n'
                    #aligned_data.append(aligned_line)
                    adf.write(datadelimiter.join(aligned_list)+'\n')
                    linecount +=1
        
        # Append new data file to aligned_data_list
        aligned_data_file_list.append(aligned_data_path_file)
        timesteps = linecount
    return aligned_data_file_list




def checktype(file_list, file_type):
    # Returns true if all files in file_list are of type file_type
    file_type_check = True
    chars = len(file_type)
    for path_file in file_list:
        if path_file[-chars:] != file_type:
            file_type_check = False
            break
    return file_type_check
