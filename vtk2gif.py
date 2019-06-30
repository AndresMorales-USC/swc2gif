# -*- coding: utf-8 -*-
"""
Created on 2018/05/22

@author: Andres Morales

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


Required Argument:
vtk_path_file: A string containing the file path and file name of the
    ".vtk" file. If using a group of sequential vtk's, the input
    vtk_path_file does NOT change, but the actual file names are
    appended with '_' and a number (starting at 0).
    (Example: 'C:/Retina simulation/morphology/NewWaveform.vtk'
     Actual File Name(s):
            NewWaveform.vtk
        or: NewWaveform_0.vtk, NewWaveform_1.vtk, NewWaveform_2.vtk)


Optional Kew Word Arguments:
save_path: A string containing the file path to where the gif(s) should
    be saved. If not given, the location of the vtk(s) will be used.
    (Example: 'C:/Retina simulation/morphology')

datatitle: Name used within the vtk for relevant cell data. Defaults
    to 'vdata'.

size: A tuple containing the horizontal and vertical dimensions of the
    image frame. Defaults to (1000,1000).

showaxes: A boolean for whether the coordinate axes in the frame is
    visible: Defaults to True.

colorBounds: A tuple or list containing the minimum and maximum values
    for the color bar. Defaults to an empty tuple, which results in
    searching for the global min/max values across all values stored in
    the cell data.

time_list: A list or tuple containing the time stamps of the cell data.
    Each element is the time (milliseconds) of each frame. If given,
    the title in the gif will use these time stamps instead of the time
    stamps/frame numbers in the vtk.
    (Example: list(range(0, 25000)))

selectframes: A list containing integers and/or sublists of integers.
    The numbers are the frames to be combined into gif(s), where each
    element of selectframes is a different gif. If not given (or set
    to an empty list), each vtk split ('model_6.vtk') will be converted
    into a separate gif. Defaults to an empty list.
    (Example: [3, [1,2,3], 4, 5, 10, [9,12,17]])

fpvtk: An int of the maximum number of frames per vtk. Setting
    value to greater than zero prevents vtk2gif from having to find the
    value itself. Defaults to 0.
    
totalframes: An int of the total number of frames accross all vtks.
    Setting value to greater than zero prevents vtk2gif from having to
    find the value itself. Defaults to 0.

startpositionoffset: A list/tuple of 3 floats that is converted into a
    numpy array. This offsets the automatically set focal point and
    camera positions.
    Default is (0.0,0.0,0.0).

startradius: A float for the change from the automatically set value for
    the distance of camera from the focal point for the first animation
    frame.
    Default is 0.

stepposition: A list/tuple of 3 floats that is converted into a numpy
    array. The list sets the change in focal point and camera position
    for each frame of animation.
    Default is (0.0,0.0,0.0).

stepradius: A float for the change in distance of camera and focal point
    each frame of animation.

angleradians: A boolean for whether angles provided are in radians.
    Defaults to False.

startelevation: A float for the starting elevation angle. Defaults to 0
    degrees from the +z-axis, where positive degrees is towards the
    +y-axis.

startazimuth: A float for the starting azimuth angle. Defaults to 0
    degrees from the +z-axis, where positive degrees is counter
    clockwise around the +y-axis.

stepelevation: A float for the elevation angle change each frame.
    Defaults to 0. Will NOT increment past -90 or 90 degrees.

stepazimuth: A float for the azimuth angle change each frame.
    Defaults to 0.
    
fps: A positive integer for the frames per second used by the gif.
    Defaults to 50.

bgcolor: Background color of animation is a tuple of rgb values.
    Defaults to black (0.0,0.0,0.0)
    
workers: A positive interger dictating the quantity of workers for
    parallel processing. Defaults to the max number of cpu's, or, if
    parallel processing isn't implemented, to 1.

Prerequisite Packages: tqdm, imageio, mayavi
pip install mayavi
pip install --upgrade numpy tqdm imageio

"""

#alternate way of running: mayavi2 -x 'C:\Users\Andres\A_testchamber\vtk2gif.py'

import os
import re
import math
from mayavi import mlab
from mayavi.sources.vtk_file_reader import VTKFileReader
from mayavi.modules.surface import Surface
from tqdm import tqdm
from numpy import array
import imageio

from multiprocessing import Pool, cpu_count


# Load VTK
def load_vtk(temp_vtk_path_file, engine, datatitle, visible=True, useCBar=True):
    # Clear current scene of possible vtk
    engine.scenes[0].children[0:1] = []
    # Add file data to engine
    #print('\nOpening: '+temp_vtk_path_file)
    src = engine.open(temp_vtk_path_file)
    
    if visible:
        # Add surface module to make data visible
        engine.add_filter(Surface(), src)
    
    # Populate list of indices of the scalars, within the vtk, that have datatitle as part of its name
    datalist = []
    for dataID, dataname in enumerate(src._cell_scalars_list):
        if datatitle in dataname:
            datalist.append(dataID)
    
    if useCBar:
        # Make colorbar visible, oriented, scaled, and positioned
        module_manager = src.children[0]
        module_manager.scalar_lut_manager.show_scalar_bar = True
        module_manager.scalar_lut_manager.scalar_bar_representation.orientation = 0
        #module_manager.scalar_lut_manager.scalar_bar.orientation = 'horizontal'
        module_manager.scalar_lut_manager.scalar_bar_representation.position2 = array([0.8, 0.1]) #size
        module_manager.scalar_lut_manager.scalar_bar_representation.position = array([0.1, 0.87]) #position of bottom-left-most point
        module_manager.scalar_lut_manager.scalar_bar.maximum_number_of_colors = 21 #limits color range for consistent gif colors
        #[0.1, 0.01] for bottom alignment
    
    return src, datalist

# Render animation frames in single engine
def render_gif_para(size, bgcolor, showaxes, gif_frame_sublist, gif_file, fpvtk, vtk_path_file_list, 
                    startpositionoffset, startradiusoffset, startelevation, startazimuth,
                    vtk_name, save_path, fps, stepposition, stepradius, stepelevation, stepazimuth,
                    datatitle, time_list, usertimeChars, minV, maxV):
            # Create new scene with set window size (which defines the gif size)
            scn = mlab.figure(size=size, bgcolor=bgcolor) #bgcolor must be tuple (0.0,0.0,0.0) - (1.0,1.0,1.0)
            # Enable view of x,y,z-axes labels
            scn.scene.show_axes = showaxes
            # Set view to +z axes
            scn.scene.z_plus_view()
            # Set clipping range
            #scn.scene.camera.clipping_range = [100,100]
            # Get engine for adding file data
            engine = mlab.get_engine()
            
            # GIF creation: Loop through the list of sublists (each sublist contains the frames for a single gif)
            # Initialize before looping
            vtknum = int(gif_frame_sublist[0]/fpvtk)
            temp_vtk_path_file = vtk_path_file_list[vtknum]
            
            # Load vtk
            src, datalist = load_vtk(temp_vtk_path_file, engine, datatitle, visible=True, useCBar=True)
            module_manager = src.children[0]
            
            # Get the auto set view angle data
            start_focal_point = scn.scene.camera.focal_point
            cam_pos = scn.scene.camera.position
            rel_cam_pos = cam_pos-start_focal_point
            start_rel_radius = (rel_cam_pos[0]**2+rel_cam_pos[1]**2+rel_cam_pos[2]**2)**0.5
        
            # Set focal point, camera radius, and angles to start values
            focal_point = start_focal_point + array(startpositionoffset)
            rel_radius = start_rel_radius + startradiusoffset
            elevation = startelevation
            azimuth = startazimuth
            new_rel_cam_pos = array([rel_radius*math.cos(elevation)*math.sin(azimuth), rel_radius*math.sin(elevation),rel_radius*math.cos(elevation)*math.cos(azimuth)])
            new_cam_pos = new_rel_cam_pos+focal_point
            scn.scene.camera.focal_point = focal_point
            scn.scene.camera.position = new_cam_pos
            
            gif_im_sublist = []
            # Loop through each frame in sublist
            for frameIDX in tqdm(gif_frame_sublist, desc='Generating GIF: '+gif_file):
                # Frame 0000 is at index 0 of vtk 0
                # Frame 0500, when fpvtk=500, is at index 0 of vtk 1
                vtkIDX = frameIDX%fpvtk
                
                # If the vtk has changed, load the new vtk
                temp_vtknum = int(frameIDX/fpvtk)
                if vtknum != temp_vtknum:
                    vtknum = temp_vtknum
                    temp_vtk_path_file = vtk_path_file_list[vtknum]
                    
                    # Load vtk
                    src, datalist = load_vtk(temp_vtk_path_file, engine, datatitle, visible=True, useCBar=True)
                    module_manager = src.children[0]
                
                
                # Change data to next selected frame
                src._cell_scalars_name = src._cell_scalars_list[datalist[vtkIDX]]
                
                # Set camera position and angle
                new_rel_cam_pos = array([rel_radius*math.cos(elevation)*math.sin(azimuth), rel_radius*math.sin(elevation), rel_radius*math.cos(elevation)*math.cos(azimuth)])
                new_cam_pos = new_rel_cam_pos+focal_point
                scn.scene.camera.focal_point = focal_point
                scn.scene.camera.position = new_cam_pos
                scn.scene.renderer.reset_camera_clipping_range()
                #scn.scene.camera.clipping_range = [100,100]
                #scn.scene.camera.compute_view_plane_normal()
                #scn.scene.render()
                
                # Step camera position and angle values for future frames
                focal_point += array(stepposition)
                rel_radius += stepradius
                if elevation+stepelevation >= math.radians(-90) and elevation+stepelevation <= math.radians(90):
                    elevation += stepelevation
                azimuth += stepazimuth
                
                # Update frame for viewing
                src.update_data()
                
                # Enforce colorbar title and limits
                padding = '_'
                framenum_padded, timepoint_padded = os.path.split(src._cell_scalars_name[len(datatitle):])
                framenum = framenum_padded.lstrip(padding)
                timepoint = timepoint_padded.lstrip(padding)
                if len(time_list)>0: # Use time provided by user
                    usertime = str(time_list[int(framenum)])
                    # Add padding for consistent title formatting/spacing
                    timeStepString = 'Time(ms): '+' '*(usertimeChars-len(usertime))+usertime
                else: # Use frame/time in vtk
                    if len(timepoint)>0: # Time data exists in vtk, use it
                        timeStepString = 'Time(ms): '+' '*(len(timepoint_padded)-len(timepoint))+timepoint
                    else: # Only frame number in vtk
                        timeStepString = 'Frame: '+' '*(len(framenum_padded)-len(framenum))+framenum
                module_manager.scalar_lut_manager.scalar_bar.title = 'Voltage (mV) at '+timeStepString
                module_manager.scalar_lut_manager.data_range = array([minV, maxV])
                
                # Append the figure to the gif animation
                f = mlab.gcf(engine=engine)
                f.scene._lift()
                gif_im_sublist.append(mlab.screenshot())
                    
            # Clear current vtk from scene
            #engine.scenes[0].children[0:1] = []
        
            # Close figure
            mlab.close()
            return gif_im_sublist



def vtk2gif(*arg, **kwargs):
    # Check if vtk file path was passed through arg
    if not len(arg) == 1:
        raise TypeError('vtk2gif expected 1 argument (vtk_path_file), got '+str(len(arg)))
    elif not isinstance(arg[0], str):
        raise TypeError('vtk2gif expected str for vtk_path_file, got '+str(type(arg[0])))
    elif not arg[0][-4:] == '.vtk':
        raise TypeError('vtk2gif expected vtk_path_file string to end in \'.vtk\'')
    else:
        vtk_path, vtk_name = os.path.split(arg[0])
    
    # Get list of vtks: either 1 vtk (vtk_name.vtk) or a set (vtk_name_###.vtk)
    vtk_path_file_list = []
    valid = re.compile(vtk_name[:-4]+'(_[0-9]+)?\.vtk')
    for afile in os.listdir(vtk_path):
        if valid.fullmatch(afile):
            vtk_path_file_list.append(os.path.join(vtk_path, afile))
    if len(vtk_path_file_list) == 0:
        raise FileNotFoundError('No such file: '+vtk_name+' at path: '+vtk_path)
    # Ensure that vtk file list is sorted properly
    dre = re.compile(r'(\d+)')
    vtk_path_file_list.sort(key=lambda l: [int(s) if s.isdigit() else s.lower() for s in re.split(dre, l)])
    numvtks = len(vtk_path_file_list)

    
    # Check if keys are defined in kwargs
    # Else set to defaults
    
    # Save path defaults to the location of the vtk's
    save_path = kwargs.get('save_path', vtk_path)
    if not isinstance(save_path, str):
        raise TypeError('vtk2gif expected str for save_path, got '+str(type(save_path)))
    
    # Title of data default: 'vdata'
    datatitle = kwargs.get('datatitle', 'vdata')
    if not isinstance(datatitle, str):
        raise TypeError('vtk2gif expected str for datatitle, got '+str(type(datatitle)))
    
    # Gif size default: (1000, 1000)
    size = kwargs.get('size', (1000, 1000))
    if not (isinstance(size, tuple) or isinstance(size, list)):
        raise TypeError('vtk2gif expected tuple or list for size, got '+str(type(size)))
    elif not len(size) == 2:
        raise TypeError('vtk2gif expected tuple or list with 2 elements for size, got '+str(len(size))+' elements')
    
    # Show axis marker defualt: True
    showaxes = kwargs.get('showaxes', True)
    if not isinstance(showaxes, bool):
        raise TypeError('vtk2gif expected bool for showaxes, got '+str(type(showaxes)))
    
    # Default frames to be turned into gif's: []
    # can be a list of integers, a list of sublists of integers, or a 
    #     string of a file path
    # file format: each line is numbers separated by spaces
    #     indicating the frames of a single gif
    selectframes=kwargs.get('selectframes', [])
    
    # Time data default: [] #Results in using frame/time info from vtk file
    time_list = kwargs.get('time_list', [])
    if not (isinstance(time_list, tuple) or isinstance(time_list, list)):
        raise TypeError('vtk2gif expected tuple or list for time_list, got '+str(type(time_list)))
    # Get max character length of timepoints
    usertimeChars = 0
    for t in time_list:
        usertimeChars = max(usertimeChars, len(str(t)))
    
    # Maximum frames per vtk default: 0 #Results in reading number of frames in 1st vtk
    fpvtk = kwargs.get('fpvtk', 0)
    if not isinstance(fpvtk, int):
        raise TypeError('vtk2gif expected int for fpvtk, got '+str(type(fpvtk)))
    elif fpvtk < 0:
        raise TypeError('vtk2gif expected int greater than or equal to zero for fpvtk, got '+str(fpvtk))
        
    # Total frames across all vtks default: 0 #Results in reading number of frames in last vtk
    totalframes = kwargs.get('totalframes', 0)
    if not isinstance(totalframes, int):
        raise TypeError('vtk2gif expected int for totalframes, got '+str(type(totalframes)))
    elif totalframes < 0:
        raise TypeError('vtk2gif expected int greater than or equal to zero for totalframes, got '+str(totalframes))
    
    # Bounds of colorbar default: () #Results in searching for data extremes
    colorBounds = kwargs.get('colorBounds', ())
    if not (isinstance(colorBounds, tuple) or isinstance(colorBounds, list)):
        raise TypeError('vtk2gif expected tuple or list for colorBounds, got '+str(type(colorBounds)))
    elif not (len(colorBounds) == 0 or len(colorBounds) == 2):
        raise TypeError('vtk2gif expected tuple or list with 2 elements for colorBounds, got '+str(len(colorBounds))+' elements')
    
    # Default offsets and movement of camera focal point: (0.0,0.0,0.0) and radius: 0.0
    startpositionoffset = kwargs.get('startpositionoffset', (0.0,0.0,0.0))
    stepposition = kwargs.get('stepposition', (0.0,0.0,0.0))
    if not (isinstance(startpositionoffset, tuple) or isinstance(startpositionoffset, list)):
        raise TypeError('vtk2gif expected tuple or list for startpositionoffset, got '+str(type(startpositionoffset)))
    elif not (len(startpositionoffset) == 3):
        raise TypeError('vtk2gif expected tuple or list with 3 elements for startpositionoffset, got '+str(len(startpositionoffset))+' elements')
    else:
        for i, v in enumerate(startpositionoffset):
            if not (isinstance(v, int) or isinstance(v, float)):
                raise TypeError('vtk2gif expected int or float for elements of startpositionoffset, got '+str(type(v))+' for element '+str(i))
    if not (isinstance(stepposition, tuple) or isinstance(stepposition, list)):
        raise TypeError('vtk2gif expected tuple or list for stepposition, got '+str(type(stepposition)))
    elif not (len(stepposition) == 3):
        raise TypeError('vtk2gif expected tuple or list with 3 elements for stepposition, got '+str(len(stepposition))+' elements')
    else:
        for i, v in enumerate(stepposition):
            if not (isinstance(v, int) or isinstance(v, float)):
                raise TypeError('vtk2gif expected int or float for elements of stepposition, got '+str(type(v))+' for element '+str(i))
    
    # Default offsets and movement of camera radius: 0.0
    startradiusoffset = kwargs.get('startradiusoffset', 0.0)
    stepradius = kwargs.get('stepradius', 0.0)
    if not (isinstance(startradiusoffset, int) or isinstance(startradiusoffset, float)):
        raise TypeError('vtk2gif expected int or float for startradiusoffset, got '+str(type(startradiusoffset)))
    if not (isinstance(stepradius, int) or isinstance(stepradius, float)):
        raise TypeError('vtk2gif expected int or float for stepradius, got '+str(type(stepradius)))
    
    # Default unit of camera angle: degrees
    angleradians =  kwargs.get('angleradians', False)
    if not isinstance(angleradians, bool):
        raise TypeError('vtk2gif expected bool for angleradians, got '+str(type(angleradians)))
    
    # Default camera elevation and azimuth angles: 0 (degrees)
    if angleradians:
        startelevation = kwargs.get('startelevation', 0)
        startazimuth = kwargs.get('startazimuth', 0)
        stepelevation = kwargs.get('stepelevation', 0)
        stepazimuth = kwargs.get('stepazimuth', 0)
    else:
        startelevation = math.radians(kwargs.get('startelevation', 0))
        startazimuth = math.radians(kwargs.get('startazimuth', 0))
        stepelevation = math.radians(kwargs.get('stepelevation', 0))
        stepazimuth = math.radians(kwargs.get('stepazimuth', 0))
    if not (isinstance(startelevation, int) or isinstance(startelevation, float)):
        raise TypeError('vtk2gif expected int or float for startelevation, got '+str(type(startelevation)))
    if not (isinstance(startazimuth, int) or isinstance(startazimuth, float)):
        raise TypeError('vtk2gif expected int or float for startazimuth, got '+str(type(startazimuth)))
    if not (isinstance(stepelevation, int) or isinstance(stepelevation, float)):
        raise TypeError('vtk2gif expected int or float for stepelevation, got '+str(type(stepelevation)))
    if not (isinstance(stepazimuth, int) or isinstance(stepazimuth, float)):
        raise TypeError('vtk2gif expected int or float for stepazimuth, got '+str(type(stepazimuth)))
    
    # Default is 50 frames per seconds of animation
    fps = kwargs.get('fps', 50)
    if not (isinstance(fps, int) or isinstance(fps, float)):
        raise TypeError('vtk2gif expected int or float for fps, got '+str(type(fps)))

    # Background color is default to black (r=0.0, g=0.0, b=0.0)
    bgcolor = kwargs.get('bgcolor', (0.0,0.0,0.0))
    if not (isinstance(bgcolor, tuple) or isinstance(bgcolor, list)):
        raise TypeError('vtk2gif expected tuple or list for bgcolor, got '+str(type(bgcolor)))
    elif not (len(bgcolor) == 3):
        raise TypeError('vtk2gif expected tuple or list with 3 elements for bgcolor, got '+str(len(bgcolor))+' elements')
    else:
        for i, v in enumerate(bgcolor):
            if not (isinstance(v, int) or isinstance(v, float)):
                raise TypeError('vtk2gif expected int or float for elements of bgcolor, got '+str(type(v))+' for element '+i)
    
    # Parse and error check selection of frames
    select_frame_list = []
    # Should be a list
    if not isinstance(selectframes, list):
        raise TypeError('vtk2gif expected a list (containing ints and/or sublists of ints) for selectframes, got '+str(type(selectframes)))
    else:
        # Check if list contains ints and/or sublists of ints, and if valid append elements to select_frame_list
        for idx, element in enumerate(selectframes):
            # Check if list conatins ints and or lists
            if not isinstance(element, int) and not isinstance(element, list):
                raise TypeError('vtk2gif expected a ints or sublists of ints for elements of selectframes, got '+str(type(selectframes[idx]))+' at element '+str(idx))
            elif isinstance(element, list):
                # Check if sublists contain ints
                for subidx, subelement in enumerate(element):
                    if not isinstance(subelement, int):
                        raise TypeError('vtk2gif expected ints for elements of sublists within selectframes, got '+str(type(selectframes[idx][subidx]))+' at ['+str(idx)+']['+str(subidx))
                select_frame_list.append(element)
            else: # Element is an int and is appended as sublist to select_frame_list
                select_frame_list.append([element])
    
    # Number of workers for parallel processing is default to the max number of cpu's
    # (unless parallel isn't implemented then it defaults to 1)
    try:
        tempworkers = cpu_count()
    except NotImplementedError:
        tempworkers = 1
    workers = kwargs.get('workers', tempworkers)
    if not isinstance(workers, int):
        raise TypeError('vtk2gif expected an int for workers, got '+str(type(workers)))
    else:
        if not workers>=1:
            raise TypeError('vtk2gif expected a positive int for workers, got '+str(workers))




    
    if fpvtk == 0 or totalframes == 0 or len(colorBounds) == 0:
        # Create new scene with set window size (which defines the gif size)
        scn = mlab.figure(size=size, bgcolor=bgcolor) #bgcolor must be tuple (0.0,0.0,0.0) - (1.0,1.0,1.0)
        # Enable view of x,y,z-axes labels
        scn.scene.show_axes = showaxes
        # Set view to +z axes
        scn.scene.z_plus_view()
        # Set clipping range
        #scn.scene.camera.clipping_range = [100,100]
        # Get engine for adding file data
        engine = mlab.get_engine()
        
        # Find fpvtk by user input (already set >0) or counting the number of frames in the 1st vtk of vtk_path_file_list
        if fpvtk == 0:
            print('Finding frames per vtk')
            # Load vtk
            src, datalist = load_vtk(vtk_path_file_list[0], engine, datatitle, visible=False, useCBar=False)
            fpvtk = len(datalist)
        
        # Find total frames in all vtk(s) if not already set by user
        if totalframes == 0:
            print('Finding total frames in vtk set')
            if numvtks == 1:
                totalframes = fpvtk
            elif numvtks > 1:
                # Load vtk
                src, datalist = load_vtk(vtk_path_file_list[-1], engine, datatitle, visible=False, useCBar=False)
                fplastvtk = len(datalist)
                totalframes = fpvtk*(numvtks-1)+fplastvtk
        
        # Set minimum and maximum values of the colorbar
        if len(colorBounds) == 0:
            print('Finding minimum and maximum data values in vtk set')
            # Find global extremes of data for colorbar bounds
            # (otherwise the bounds would be determined by the local min/max
            # values and change with each frame)
            
            # Initialize min/max
            minV = float('Inf')
            maxV = -float('Inf')
            
            # Loop through each vtk file seaching for data extremes
            for temp_vtk_path_file in tqdm(vtk_path_file_list, desc='Finding data minimum and maximum'):
                # Load vtk
                src, datalist = load_vtk(temp_vtk_path_file, engine, datatitle, visible=False, useCBar=False)
                module_manager = src.children[0]
                
                # Get global min/max values of this vtk
                # Loop through all frames
                for vdataID in tqdm(datalist, desc='Finding local min/max of '+os.path.split(temp_vtk_path_file)[1]):
                    # Change data to next step
                    src._cell_scalars_name = src._cell_scalars_list[vdataID]
                    # Update frame
                    src.update_data()
                    # Get min/max values of this frame
                    tempMin = module_manager.scalar_lut_manager.data_range[0]
                    tempMax = module_manager.scalar_lut_manager.data_range[1]
                    if tempMin < minV:
                        minV = tempMin
                    if tempMax > maxV:
                        maxV = tempMax
            print('\nData Extremes:')
            print([minV, maxV])
        
        # Close figure
        mlab.close()
    
        
    # If no frames are selected, select frames such that each vtk is converted into a gif
    if len(select_frame_list) == 0:
        for vtknum in range(numvtks-1):
            select_frame_list.append(list(range(vtknum*fpvtk,(vtknum+1)*fpvtk)))
        select_frame_list.append(list(range((numvtks-1)*fpvtk, totalframes)))
    
    # Check if values are within the bounds of the vtk(s)
    if min(min(select_frame_list)) < 0 or max(max(select_frame_list)) > totalframes:
        raise TypeError('vtk2gif expected selectframes to be within the bounds of vtk(''s)')
        
    if len(colorBounds) > 0:
        minV, maxV = colorBounds
    
    
    
    
    pool = Pool(processes=workers)
    # Loop through frame sets each to be converted to a single gif
    for gif_frame_list in tqdm(select_frame_list, desc='Generating GIF(s)'):
        # Split up selected frames (of singular gif) among multiple workers
        gif_frame_sublists = [gif_frame_list[i:i + workers] for i in range(0, len(gif_frame_list), workers)]
        # Establish starting positions for each worker
        startpositionoffset_list = []
        startradiusoffset_list = []
        startelevation_list = []
        startazimuth_list = []
        steps = 0
        for w_idx in range(workers):
            startpositionoffset_list.append(tuple([startpositionoffset[i]+stepposition[i]*steps for i in range(3)]))
            startradiusoffset_list.append(startradiusoffset+stepradius*steps)
            startelevation_list.append(startelevation+stepelevation*steps)
            startazimuth_list.append(startazimuth+stepazimuth*steps)
            steps = steps + len(gif_frame_sublists[w_idx])
        
        # Generate gif name
        gif_file = vtk_name[:-4]+' - frames '+str(gif_frame_list[0])+'-'+str(gif_frame_list[-1])+'.gif'
        gif_path_file = os.path.join(save_path, gif_file)
        
        # Get images from workers
        gif_im_list = pool.starmap(render_gif_para,
                                [(size, bgcolor, showaxes, gif_frame_sublists[w_idx], gif_file, fpvtk, vtk_path_file_list, 
                                startpositionoffset_list[w_idx], startradiusoffset_list[w_idx],
                                startelevation_list[w_idx], startazimuth_list[w_idx],
                                vtk_name, save_path, fps, stepposition, stepradius, stepelevation, stepazimuth,
                                datatitle, time_list, usertimeChars, minV, maxV)
                                for w_idx in range(workers)])
        
        
        # Export frames
        imageio.mimwrite(gif_path_file, [image for sublist in gif_im_list for image in sublist], fps=fps, subrectangles=False) #, mode='I'
        
        # with imageio.get_writer(gif_path_file, mode='I', fps=fps) as writer:
            # for gif_im_list in gif_im_listlist:
                # for im in gif_im_list:
                    # writer.append_data(im)                 
        