import os
from swc2gif.swc2vtk import swc2vtk
from swc2gif.vtk2gif import vtk2gif
from time import time

# Can run with (ipython has cleanest output, python is fastest)
# python AD_gifs_Vm.py
# ipython AD_gifs_Vm.py
# mayavi2 -x 'AD_gifs_Vm.py'

# __main__ is required for parallel processing not to get stuck in a loop
if __name__ == '__main__':
    
    t0 = time()
    
    ## Set file paths
    save_path = 'C:/Users/Andres/A_testchamber/working data'
    
    # Create list of all swc morphology files in folder swc_path
    swc_path = 'C:/Users/Andres/A_testchamber/Retina simulation/morphology'
    morphology_list = []
    for afile in os.listdir(swc_path):
        if afile.endswith(".swc"):
            swc_path_file = os.path.join(swc_path, afile)
            morphology_list.append(swc_path_file)

    voltage_path_file = 'C:/Users/Andres/A_testchamber/simulation data/sectionValues_inc1.v'
    datadelimiter =  ' '
    coords_path_file = 'C:/Users/Andres/A_testchamber/simulation data/sectionCoords_inc1.txt'
    coordsdelimiter = '	'
    time_path_file = 'C:/Users/Andres/A_testchamber/simulation data/sectionTimes_inc1.txt'

    ## Set swc2vtk() parameters
    vtk_name='NeuroModel.vtk'
    datatitle ='vdata'
    spherediv=8
    cyldiv=5
    fpvtk = 500 #Frames per vtk

    # Convert swc(s) into vtk(s)
    minV, maxV, totalframes = swc2vtk(morphology_list, voltage_path_file, coords_path_file,
                                   datadelimiter=datadelimiter, coordsdelimiter=coordsdelimiter,
                                   time_path_file=time_path_file,
                                   save_path=save_path, vtk_name=vtk_name, datatitle=datatitle,
                                   spherediv=spherediv, cyldiv=cyldiv, fpvtk=fpvtk)

    t1 = time()
    print('\nSWC2VTK Time (s): ')
    print(t1-t0)
    print('\n')

    # Set vtk2gif() parameters
    vtk_path = os.path.join(save_path, vtk_name)
    # Select the frames to be turned into gif(s)
    # To select 100-350 ms (and data is 0.2 ms/frame), we would need to select frames 500-1750
    selectframes = [list(range(505, 600))] #example alternative (2 gifs each 1 frame): [[200],[201]]
    colorBounds = (minV, maxV)
    fps=10 # Playback speed (frames per second)
    startpositionoffset = [0,50,0]
    startradiusoffset = -200
    
    # Optional electrode arguments
    elecX = [70-150, -150]
    elecY = [70+150, 150]
    elecZ = [10, 0]
    elecradius = 1
    eleccolor = (1,1,1)
    
    # # Other optional animation keyword arguments
    # stepposition = [0,-50/1250, 0]
    # stepradius = 200/1250
    # startelevation=50
    # stepelevation=-0.15
    # startazimuth=130
    # stepazimuth=0.3

    # Convert vtk(s) and time data into animated gif(s)
    
    # Convert vtk(s) into gif(s)
    vtk2gif(vtk_path, selectframes=selectframes, datatitle=datatitle, save_path=save_path,
            colorBounds=colorBounds, fps=fps, fpvtk=fpvtk, totalframes=totalframes,
            startpositionoffset=startpositionoffset, startradiusoffset=startradiusoffset,
            elecX=elecX, elecY=elecY, elecZ=elecZ, elecradius=elecradius, eleccolor=eleccolor)
            #, stepposition=stepposition, stepradius=stepradius, startelevation=startelevation, stepelevation=stepelevation, startazimuth=startazimuth, stepazimuth=stepazimuth)

    t2 = time()
    print("\nVTK2GIF Time (s): ")
    print(t2-t1)
    print('\n')
