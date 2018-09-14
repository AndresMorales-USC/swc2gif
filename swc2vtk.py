# -*- coding: utf-8 -*-
"""
Created on 2018/05/22

@author: Andres Morales

main(swc_file_path(s), [data_file_path(s)], [coordinate_file_path(s)],
    [save_path=swc_path, vtkname='model.vtk', datatitle='vdata',
    datadelimiter=' ', coordsdelimiter='	', maxframes=float('Inf')])

returns number of vtk's saved and the min and max data values of all
    frames of all vtk's


Required Argument:
swc_file_path(s): A string containing the file path and file name,
    or tuple/list of strings containing the file paths and file names,
    of the ".swc" file(s).
    (Example 0: 'C:\\Retina simulation\\morphology\\606.swc'
     Example 1: ('C:\\Retina simulation\\morphology\\606.swc',
                'C:\\Retina simulation\\morphology\\1010.swc')
            or: ['C:\\Retina simulation\\morphology\\606.swc',
                'C:\\Retina simulation\\morphology\\1010.swc'])


Optional Arguments:
data_file_path(s): A string containing the file path and file name,
    or tuple/list of strings containing the file paths and file names,
    of the ".txt" data file(s). Each data file is information to be
    mapped to a specific swc (swc & data file lists must be of equal
    size and ordering). Or, data may be in a single file that requires
    coordinate information to map onto the morphologies.
    (Example 1a: ['C:\\Retina simulation\\morphology\\606_data.txt',
                  'C:\\Retina simulation\\morphology\\1010_data.txt']
     Example 1b: 'C:\\Retina simulation\\morphology\\voltage_data.txt')

coordinate_file_path(s): A string containing the file path and file
    name of the ".txt" coordinate file.
    (Example 1b: 'C:\\Retina simulation\\morphology\\coords_data.txt')


Optional Kew Word Arguments:
save_path: A string containing the file path of where to save generated
    files (alinged data files or ".vtk" files). Save_path defaults to 
    the same folder as where the ".swc" files are read.
    (Example: 'C:\\Retina simulation\\outputs')

vtkname: A string containing the file name of the ".vtk" file to be
    generated. vtkname, if not given, is set to 'model.vtk'.

datatitle: A string containing the name to label data within the ".vtk"
    file. It defaults to 'vdata'.

datadelimiter: A string used to split input from the data file. The
    default data delimiter is ' '.

coordsdelimiter: A string used to split input from the coordinate file.
    The default coordinate data delimiter is '	'.

maxframes=float('Inf'): An integer for the maximum number of frames to
    be written to a single ".vtk". If the number of frames in the data
    file exceedes maxframes, multiple ".vtk"s will be written. Each vtk
    containing at most the frame quota and sequentially labeled with
    the prefix '_' and a number (starting at 0). The default value is
    actually a float for infinity.

"""

import gifgen.vtkgen as vtkgen
import os

def main(*args, **kwargs):
    # Check if swc, data, and coordinate file paths were passed through arg
    # If none ask for all paths
    # (data propt can be canceled, vtk just won't have data)
    # (coords can be canceled if data is already aligned)
    # If only 1, assume data and coords were intentionally not included
    # If only 2, assume data was already aligned
    # If more than 3 will ignore extra args
    if len(args) == 0:
        print('ERROR: Not enough arguments [minimum: list of swc''s]')
    else: # SWC arguments are contained in args
        # Place file path strings into a list from whatever format it was submitted
        if isinstance(args[0], basestring):
            swc_list = [args[0]]
        elif isinstance(args[0], tuple):
            swc_list = list(args[0])
        elif isinstance(args[0], list):
            swc_list = args[0]
        else:
            print('Warning: Unknown input for swc file(s)')
        
        data_file_list = []
        if len(args) > 1: # Data arguments are contained in args
            if isinstance(args[1], basestring):
                data_file_list = [args[1]]
            elif isinstance(args[1], tuple):
                data_file_list = list(args[1])
            elif isinstance(args[1], list):
                data_file_list = args[1]
            else:
                print('Warning: Unknown input for data file(s)')
            
        coords_file_path = ''
        if len(args) > 2: # Coordinate argument is contained in args
            if isinstance(args[2], basestring):
                coords_file_path = args[2]
            else:
                print('Warning: Unknown input for coordinate file(s)')
            
        if len(args) > 3: # Extra arguments are contained in args
            print('Warning: Extra arguments recieved beyond SWC, Data, and Coordinate files')
    
    
    # Check if number of files for each type is consistent
    # number of data files can be 1 with an associated coordinate file
    # or one data file for each swc and no coordinate files
    
    # If coordinate TXT's are used, check if number is same as data V's
    if len(coords_file_path) > 0 and len(data_file_list) != 1:
        print('Warning: Mismatch in number of data V and coordinate TXT files')
    
    # If data V's are used, check if number is appropriate for SWC's
    if len(data_file_list)>1 and len(data_file_list) != len(swc_list):
        print('Warning: Mismatch in number of SWC and data V files')
    
    
    # Check file types for all args
    # Check if all SWC's have swc as file types
    if not checktype(swc_list, 'swc'):
        print('Warning: One or more SWC\'s do NOT have a \'.swc\' file type')
        
    # Check if all data V's have v as file types
    if not checktype(data_file_list, 'v'):
        print('Warning: One or more data V\'s do NOT have a \'.v\' file type')
    
    # Check if all coordinate TXT's have txt as file types
    if not checktype([coords_file_path], 'txt'):
        print('Warning: One or more coordinate TXT\'s do NOT have a \'.txt\' file type')
    
    
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
    
    
    # Initialize VTK, add and convert SWC file(s)
    head, tail = os.path.split(swc_list[0])
    vtk_file_path = head+'\\'+vtk_name
    vtk = vtkgen.VtkGenerator()
    for swc_file_path in swc_list:
        vtk.add_swc(swc_file_path)
    #vtk.convert_swc()
    
    # Check if data needs to be converted
    if len(coords_file_path) > 0:
        # Convert all data using SWC's and coordinate TXT's
        from scipy import spatial
        from tqdm import tqdm
        
        # Extract data coordinates from coords file
        tempHead, coords_file = os.path.split(coords_file_path)
        data_coords = []
        with open(coords_file_path, 'r') as f:
            # Read each line of the file as a string, split into a list of strings, convert to a
            #     list of numbers, and finally append single postion list to list of all positions
            for line in f:
                str_list = line.rstrip().split(coordsdelimiter)
                float_list = []
                for val in str_list:
                    try:
                        float_list.append(float(val))
                    except:
                        print(val)
                data_coords.append(float_list)
        
        # Generate KDTree
        tree = spatial.cKDTree(data_coords)
        
        data_file_path = data_file_list[0] # If coordinates need conversion then only one data file should be in list
        
        aligned_data_file_list = []
        timesteps=0
        # Loop through SWC files for conversion
        for file_idx, swc_file_path in enumerate(tqdm(swc_list, desc='Aligning Data to SWC\'s')):
            tempHead, swc_file = os.path.split(swc_file_path)
            
            # Iterate through each swc compartment and generate reference list
            compartment_coords_idx = []
            for swc_compartment in vtk.swc_list[file_idx].data.values():
                # Find index of nearest neighbor in coordinates and append to list
                distance, idx = tree.query(swc_compartment['pos'])
                compartment_coords_idx.append(idx)
            
            # Generate aligned data
            # If no save path was specified use the swc's path
            if len(save_path) == 0:
                aligned_data_file_path = swc_file_path[:-4]+'_aligned_data.v'
            else:
                aligned_data_file_path = os.path.join(save_path, swc_file[:-4]+'_aligned_data.v')
            with open(aligned_data_file_path, 'w') as adf:
                # Iterate through data file where each line is a different time step
                with open(data_file_path, 'r') as df:
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
            aligned_data_file_list.append(aligned_data_file_path)
            timesteps = linecount
    else: # Data is already aligned
        aligned_data_file_list = data_file_list

    # Check if data needs to be added
    if len(aligned_data_file_list) > 0:
        # Add all data to VTK and write it to the same folder as the first SWC
        for data_file_path in aligned_data_file_list:
            vtk.add_timeddatafile(data_file_path)
        numvtks, minV, maxV = vtk.write_vtk(vtk_file_path, datatitle=datatitle, delimiter=datadelimiter, maxframes=maxframes)
        
    else: # No data to add, just write the VTK to the same folder as the first SWC
        numvtks, minV, maxV = vtk.write_vtk(vtk_file_path)
    
    print('\nNumber of VTK''s: '+str(numvtks))
    print('Data Bounds [Min, Max]: ['+str(minV)+', '+str(maxV)+']')
    vtksplitdata_path = head+'\\'+vtk_name[:-4]+'_splitdata.txt'
    with open(vtksplitdata_path, 'w') as f:
        f.write('Number of VTK''s: '+str(numvtks)+'\n')
        f.write('Data Bounds [Min, Max]: ['+str(minV)+', '+str(maxV)+']')
    return [numvtks, minV, maxV]





def checktype(file_list, file_type):
    # Returns true if all files in file_list are of type file_type
    file_type_check = True
    chars = len(file_type)
    for file_path in file_list:
        if file_path[-chars:] != file_type:
            file_type_check = False
            break
    return file_type_check

if __name__ == '__main__':
    # Script is being run directly (as opposed to imported)
    # Request input of swc list, v-data, and coordinate data
    import Tkinter, tkFileDialog
    import sys
    
    root = Tkinter.Tk()
    root.withdraw()

    print('\nSelect SWC file(s)')
    print('Cancel after all SWC\'s have been opened')
    sys.stdout.flush()
    swc_list = []
    loopExit = False
    while not loopExit:
        swc_file_path = tkFileDialog.askopenfilename()
        if swc_file_path == '':
            loopExit = True
        else:
            swc_list.append(swc_file_path)
            print('Opened: ' + swc_file_path)
            sys.stdout.flush()

    print('\nSelect data V file(s)')
    print('Cancel after all data V\'s have been opened')
    print('Or cancel before selecting first data V to include no data')
    sys.stdout.flush()
    data_file_list = []
    loopExit = False
    while not loopExit:
        data_file_path = tkFileDialog.askopenfilename()
        if data_file_path == '':
            loopExit = True
        else:
            data_file_list.append(data_file_path)
            print('Opened: ' + data_file_path)
            sys.stdout.flush()
    
    print('\nSelect coordinate TXT file')
    print('Select coordinate TXT to be opened')
    print('Or cancel before selecting first coordinate TXT if data is already aligned')
    sys.stdout.flush()
    coords_file_path = tkFileDialog.askopenfilename()
    main(swc_list, data_file_list, coords_file_path)
