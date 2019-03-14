# -*- coding: utf-8 -*-
"""
Created on 2018/05/22

@author: Andres Morales

vtk2gif(vtk_path_file,
    [
    save_path=vtk_path, datatitle='vdata', size=(1000,1000),
    showaxes=True, colorBounds=(), time_path_file='', numvtks=1,
    selectframes=[] or '', maxframes=float('Inf'), angleradians=False,
    startelevation=0, startazimuth=0, stepelevation=0, stepazimuth=0,
    bgcolor=(0.0,0.0,0.0)
    ])

Takes VTK file(s) containing vdata, where each vdata## is a different
    frame of the produced GIF(s).


Required Argument:
vtk_file_path: A string containing the file path and file name of the
    ".vtk" file. If using a group of sequential vtk's, the
    vtk_file_path does NOT change, but the actual file names are
    appended with '_' and a number (starting at 0).
    (Example: 'C:\\Retina simulation\\morphology\\NewWaveform.vtk'
     Actual File Name(s):
            NewWaveform.vtk
        or: NewWaveform_0.vtk, NewWaveform_1.vtk, NewWaveform_2.vtk)


Optional Kew Word Arguments:
save_path: A string containing the file path to where the gif(s) should
    be saved. If not given, the location of the vtk(s) will be used.
    (Example: 'C:\\Retina simulation\\morphology')

datatitle: Name used within the vtk for relevant cell data. Defaults
    to 'vdata'.

size: A tuple containing the horizontal and vertical dimensions of the
    image frame. Defaults to (1000,1000).

showaxes: A boolean for whether the coordinate axes in the frame is
    visible: Defaults to True.

colorBounds: A tuple containing the minimum and maximum values for the
    color bar. Defaults to an empty tuple, which results in searching
    for the global min/max values across all values stored in the cell
    data.

time_path_file: A string containing the file path and name of the '.txt'
    document with the time stamps of the data. Each line of the
    document is the time (seconds) of each frame. If given, the title
    in the gif will use the time stamp instead of the frame number.
    (Example: 'C:\\Retina simulation\\morphology\\sectionTimes.txt')

numvtks: An integer for when using multiple sequential vtk's. If set to
    greater than 1, vtk file names are expected to end with '_#.vtk'.
    Otherwise, the vtk file name is expected to be the same as provided
    in the vtk_file_path. Defaults to 1.

selectframes: Can be a string of a file path, a list of integers, or a
    list of sublists of integers. The numbers are the frames to be
    combined into gif(s), where each sublist is a different gif. If a
    string, the file format expected is that the lines contain integers
    separated by spaces, and each line represents a sublist. If not
    given (or set to an empty list/string), each vtk will be converted
    into a separate gif. Defaults to an empty list.
    (Example: 'C:\\Retina simulation\\morphology\\selectFrames.txt'
          or: [3,4,5,10,12]
          or: [[1,2,3], [1,3,5], [9,12,17]])

maxframes: A float of the maximum number of frames in each vtk. MUST be
    set to an positive non-zero integer, if more than 1 vtk is used.
    Defaults to infinity.

startpositionoffset: A list of 3 floats that is converted into a numpy
    array. This offsets the automatically set focal point and camera
    positions.
    Default is [0.0,0.0,0.0].

startradius: A float for the change from the automatically set value for
    the distance of camera from the focal point for the first animation
    frame.
    Default is 0.

stepposition: A list of 3 floats that is converted into a numpy array.
    The list sets the change in focal point and camera position for
    each frame of animation.
    Default is [0.0,0.0,0.0].

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

Prerequisite Packages: tqdm, imageio, mayavi
pip install mayavi
pip install --upgrade numpy tqdm imageio

"""

#alternate way of running: mayavi2 -x 'C:\Users\Andres\A_testchamber\vtk2gif.py'

import os
import math
from mayavi import mlab
from mayavi.sources.vtk_file_reader import VTKFileReader
from mayavi.modules.surface import Surface
from tqdm import tqdm
from numpy import array
import imageio

def vtk2gif(*arg, **kwargs):
    
    # Check if vtk file path was passed through arg
    if len(arg) == 0:
        print('ERROR: Not enough arguments [minimum: string of vtk file path]')
    else:
        vtk_path_file = arg[0]
        vtk_path, vtk_name = os.path.split(vtk_path_file)
    
    
    # Check if keys are defined in kwargs
    # Else set to defaults
    
    # Save path defaults to the location of the vtk's
    save_path = kwargs.get('save_path', vtk_path)
    
    # Title of data default: 'vdata'
    datatitle = kwargs.get('datatitle', 'vdata')
    
    # Gif size default: (1000, 1000)
    size = kwargs.get('size', (1000, 1000))
    
    # Show axis marker defualt: True
    showaxes = kwargs.get('showaxes', True)
    
    # Bounds of colorbar default: ()
    colorBounds = kwargs.get('colorBounds', ())
    
    # Time data path default: ''
    time_path_file = kwargs.get('time_path_file', '')
    
    # Default number of vtk's to be visualized: 1
    numvtks = kwargs.get('numvtks', 1)
    
    # Default frames to be turned into gif's: []
    # can be a list of integers, a list of sublists of integers, or a 
    #     string of a file path
    # file format: each line is numbers separated by spaces
    #     indicating the frames of a single gif
    selectframes=kwargs.get('selectframes', [])
    
    # Maximum frames per vtk default: infinite
    maxframes = kwargs.get('maxframes', float('Inf'))
    
    # Default camera focal point and radius offsets and movement: 0
    startpositionoffset = array(kwargs.get('startpositionoffset', [0.0,0.0,0.0]))
    startradiusoffset = kwargs.get('startradiusoffset', 0)
    stepposition = array(kwargs.get('stepposition', [0.0,0.0,0.0]))
    stepradius = kwargs.get('stepradius', 0)
    
    # Default unit of camera angle: degrees
    angleradians =  kwargs.get('angleradians', False)
    
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
    
    # Default is 50 frames per seconds of animation
    fps = kwargs.get('fps', 50)

    # Background color is default to black (r=0.0, g=0.0, b=0.0)
    bgcolor = kwargs.get('bgcolor', (0.0,0.0,0.0))



    # If time data exists, load into an array
    time_data = []
    useTime = len(time_path_file) > 0
    if useTime:
        timeChars = 0
        with open(time_path_file, 'r') as f:
            for line in f:
                time_data.append(line.rstrip())
                temp_len = len(line)
                if temp_len > timeChars:
                    timeChars = temp_len
                
    
    # Parse and error check selection of frames
    select_frame_list = []
    if len(selectframes)>0:
        # If string pointing to file, load into a list of sublists of ints
        if isinstance(selectframes, str):
            with open(selectframes, 'r') as f:
                for line in f:
                    string_list = line.rstrip().split(' ')
                    int_list = []
                    for frame in string_list:
                        int_list.append(int(frame))
                    select_frame_list.append(int_list)
        # Else is a list
        elif isinstance(selectframes, list):
            # Check if list contains ints, and append to selec_frame_list
            if isinstance(selectframes[0], int):
                select_frame_list.append(selectframes)
            # Check if list conatins sublists of ints
            elif isinstance(selectframes[0], list):
                if isinstance(selectframes[0][0], int):
                    select_frame_list = selectframes
                else:
                    print('ERROR: selectframes is list of sublists, but doesn''t contain ints')
            else:
                print('ERROR: selectframes is a list, but doesn''t contain ints or sublists containing ints')
        else:
            print('ERROR: selectframes must be a string pointing to a file, list of ints or, a list of sublists containing ints')
        
        if min(min(select_frame_list)) < 1 or max(max(select_frame_list)) > maxframes*numvtks:
            print('ERROR: one or more frames selected are out of bounds of vtk(''s)')
    
    if (numvtks>1) and (maxframes==float('Inf')):
        print('ERROR: more than 1 vtk, but frames per vtk (maxframes) is infinite')
    
    
    
    
    # Load VTK
    def _load_vtk(*load_args, **useCBar):
        temp_vtk_path_file, engine = load_args
        useCBar = kwargs.get('useCBar', True)
        # Clear current scene of possible vtk
        engine.scenes[0].children[0:1] = []
        # Add file data to engine
        src = engine.open(temp_vtk_path_file)
        # Add surface module to make data visible
        engine.add_filter(Surface(), src)
        
        # Populate list of indices of the scalars, within the vtk, that have datatitle as part of its name
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
        
        return src
    
    
    
    
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

    
    
    # To delete current vtk from the scene
    #scene = engine.scenes[0]
    #scene.children[0:1] = []
    
    # Set minimum and maximum values of the colorbar
    if len(colorBounds) > 0:
        minV, maxV = colorBounds
    else:
        # Find global extremes of data for colorbar bounds
        # (otherwise the bounds would be determined by the local min/max
        # values and change with each frame)
        
        # Initialize min/max
        minV = float('Inf')
        maxV = -float('Inf')
        
        
        # Loop through each vtk file seaching for data extremes
        for vtknum in tqdm(list(range(numvtks)), desc='Finding data minimum and maximum'):
            # Set vtk name in case of many vtks
            if numvtks == 1:
                temp_vtk_file = vtk_name
                temp_vtk_path_file = vtk_path_file
            else:
                temp_vtk_file = vtk_name[:-4]+'_'+str(vtknum)+vtk_name[-4:]
                temp_vtk_path_file = os.path.join(vtk_path, temp_vtk_file)
            
            # Load vtk
            datalist = []
            src = _load_vtk(temp_vtk_path_file, engine, useCBar=False)
            module_manager = src.children[0]
            
            # Get global min/max values of this vtk
            # Loop through all frames
            for vdataID in tqdm(datalist, desc='Finding local min/max of '+temp_vtk_file):
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
    
    
    
    
    # If no frames are selected, convert each vtk into a gif
    if len(select_frame_list) == 0:
        for vtknum in tqdm(list(range(numvtks)), desc='Generating GIF(s)'):
            # Set vtk name in case of one or many vtks
            if numvtks == 1:
                temp_vtk_file = vtk_name
                temp_vtk_path_file = vtk_path_file
            else:
                temp_vtk_file = vtk_name[:-4]+'_'+str(vtknum)+vtk_name[-4:]
                temp_vtk_path_file = os.path.join(vtk_path, temp_vtk_file)
            gif_file = temp_vtk_file[:-3]+'gif'
            gif_path_file = os.path.join(save_path, gif_file)
            
            # Load vtk
            datalist = []
            src = _load_vtk(temp_vtk_path_file, engine, useCBar=True)
            module_manager = src.children[0]
            
            # Calculate the number of characters used for the number of the last frame (for future typography)
            frameChars = sum(c.isdigit() for c in src._cell_scalars_list[datalist[0]])
            
            # Get the auto set view angle data
            start_focal_point = scn.scene.camera.focal_point
            cam_pos = scn.scene.camera.position
            rel_cam_pos = cam_pos-start_focal_point
            start_rel_radius = (rel_cam_pos[0]**2+rel_cam_pos[1]**2+rel_cam_pos[2]**2)**0.5
            
            # Set focal point, camera radius, and angles to start values
            focal_point = start_focal_point + startpositionoffset
            rel_radius = start_rel_radius + startradiusoffset
            elevation = startelevation
            azimuth = startazimuth
            
            # Start gif creation and then iterate through datalist
            with imageio.get_writer(gif_path_file, mode='I', fps=fps) as writer:
                for frameIDX in tqdm(datalist, desc='Generating GIF: '+gif_file):
                    
                    # Change data to next vdata step
                    src._cell_scalars_name = src._cell_scalars_list[frameIDX]
                    
                    # Set camera angle
                    new_rel_cam_pos = array([rel_radius*math.cos(elevation)*math.sin(azimuth), rel_radius*math.sin(elevation), rel_radius*math.cos(elevation)*math.cos(azimuth)])
                    new_cam_pos = new_rel_cam_pos+focal_point
                    scn.scene.camera.focal_point = focal_point
                    scn.scene.camera.position = new_cam_pos
                    scn.scene.renderer.reset_camera_clipping_range()
                    #scn.scene.camera.clipping_range = [100,100]
                    #scn.scene.camera.compute_view_plane_normal()
                    #scn.scene.render()
                    
                    # Step camera angle values for future frames
                    focal_point += stepposition
                    rel_radius += stepradius
                    if elevation+stepelevation >= math.radians(-90) and elevation+stepelevation <= math.radians(90):
                        elevation += stepelevation
                    azimuth += stepazimuth
                    
                    # Update frame for viewing
                    src.update_data()
                    
                    # Enforce colorbar title and limits
                    framenum = int(src._cell_scalars_name[-frameChars:])
                    if useTime: # Use time stamps from file
                        timeIDX = framenum-1 # Frames start counting at 1, but the time index starts at 0
                        timeStep = time_data[timeIDX]
                        if len(timeStep) < timeChars: # Add padding for consistent title formatting/spacing
                            timeStepString = 'Time(ms): '+' '*(timeChars-len(timeStep))+timeStep
                        else:
                            timeStepString = 'Time(ms): '+timeStep
                    else: # Use frame number
                        framenum_str = str(framenum)
                        if len(framenum_str) < frameChars:
                            timeStepString = 'Frame: '+' '*(frameChars-len(framenum_str))+framenum_str
                        else:
                            timeStepString = 'Frame: '+str(framenum_str)
                    module_manager.scalar_lut_manager.scalar_bar.title = 'Voltage (mV) at '+timeStepString
                    module_manager.scalar_lut_manager.data_range = array([minV, maxV])
                    
                    # Append the figure to the gif animation
                    f = mlab.gcf()
                    f.scene._lift()
                    writer.append_data(mlab.screenshot())
                        
    else: # Loop through the list of sublists (each sublist contains the frames for a single gif)
        # Initialize before looping
        # Set vtk name in case of one or many vtks
        if numvtks == 1:
            temp_vtk_path_file = vtk_path_file
        else:
            vtknum = int(select_frame_list[0][0]/maxframes)
            temp_vtk_path_file = vtk_path_file[:-4]+'_'+str(vtknum)+vtk_path_file[-4:]
        
        # Load vtk
        datalist = []
        src = _load_vtk(temp_vtk_path_file, engine, useCBar=True)
        module_manager = src.children[0]
        
        # Calculate the number of characters used for the number of the last frame (for future typography)
        frameChars = sum(c.isdigit() for c in src._cell_scalars_list[datalist[0]])
        
        # Get the auto set view angle data
        start_focal_point = scn.scene.camera.focal_point
        cam_pos = scn.scene.camera.position
        rel_cam_pos = cam_pos-start_focal_point
        start_rel_radius = (rel_cam_pos[0]**2+rel_cam_pos[1]**2+rel_cam_pos[2]**2)**0.5
        
        # Loop through lists of sublists
        for gif_list in tqdm(select_frame_list, desc='Generating GIF(s)'):
            # Set focal point, camera radius, and angles to start values
            focal_point = start_focal_point + startpositionoffset
            rel_radius = start_rel_radius + startradiusoffset
            elevation = startelevation
            azimuth = startazimuth
            new_rel_cam_pos = array([rel_radius*math.cos(elevation)*math.sin(azimuth), rel_radius*math.sin(elevation),rel_radius*math.cos(elevation)*math.cos(azimuth)])
            new_cam_pos = new_rel_cam_pos+focal_point
            scn.scene.camera.focal_point = focal_point
            scn.scene.camera.position = new_cam_pos
            
            # Generate gif name and start exporting frames
            gif_file = vtk_name[:-4]+' - frames '+str(gif_list[0])+'-'+str(gif_list[-1])+'.gif'
            gif_path_file = os.path.join(save_path, gif_file)
            with imageio.get_writer(gif_path_file, mode='I', fps=fps) as writer:
                
                # Loop through each frame in sublist
                #for vdataID in tqdm(datalist[skipframes:], desc='Converting '+vtk_name+' to GIF'):
                for frameIDX in tqdm(gif_list, desc='Generating GIF: '+gif_file):
                    # If only one vtk is used, then frame selection is easy
                    if numvtks == 1:
                        vtkIDX = frameIDX
                    else:
                        # Determine which vtk contains the frame and the index of the frame within the vtk
                        # Frame 0001 is at index 0
                        if maxframes == float('Inf'):
                            vtkIDX = (frameIDX-1)
                        else:
                            vtkIDX = (frameIDX-1)%maxframes
                        
                        # Also, if the vtk has changed, load the new vtk
                        temp_vtknum = int((frameIDX-1)/maxframes)
                        if vtknum != temp_vtknum:
                            vtknum = temp_vtknum
                            temp_vtk_path_file = vtk_path_file[:-4]+'_'+str(vtknum)+vtk_path_file[-4:]
                            
                            # Load vtk
                            datalist = []
                            src = _load_vtk(temp_vtk_path_file, engine, useCBar=True)
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
                    focal_point += stepposition
                    rel_radius += stepradius
                    if elevation+stepelevation >= math.radians(-90) and elevation+stepelevation <= math.radians(90):
                        elevation += stepelevation
                    azimuth += stepazimuth
                    
                    # Update frame for viewing
                    src.update_data()
                    
                    # Enforce colorbar title and limits
                    framenum = int(src._cell_scalars_name[-frameChars:])
                    if useTime: # Use time stamps from file
                        timeIDX = framenum-1 # Frames start counting at 1, but the time index starts at 0
                        timeStep = time_data[timeIDX]
                        if len(timeStep) < timeChars: # Add padding for consistent title formatting/spacing
                            timeStepString = 'Time(ms): '+' '*(timeChars-len(timeStep))+timeStep
                        else:
                            timeStepString = 'Time(ms): '+timeStep
                    else: # Use frame number
                        framenum_str = str(framenum)
                        if len(framenum_str) < frameChars:
                            timeStepString = 'Frame: '+' '*(frameChars-len(framenum_str))+framenum_str
                        else:
                            timeStepString = 'Frame: '+str(framenum_str)
                    module_manager.scalar_lut_manager.scalar_bar.title = 'Voltage (mV) at '+timeStepString
                    module_manager.scalar_lut_manager.data_range = array([minV, maxV])
                    # Append the figure to the gif animation
                    f = mlab.gcf()
                    f.scene._lift()
                    writer.append_data(mlab.screenshot())
                    
            # Clear current vtk from scene
            #engine.scenes[0].children[0:1] = []
        
    # Close figure
    mlab.close()
