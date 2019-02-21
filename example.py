import os
from swc2gif.mapswc import mapswc
from swc2gif.swc2vtk import swc2vtk
from swc2gif.vtk2gif import vtk2gif
from time import time

# Can run with (python is fastest, ipython has cleanest output)
# python -m swc2gif.example.py
# ipython -m swc2gif.example.py
# ipython swc2gif/example.py
# mayavi2 -x 'swc2gif/example.py'

t0 = time()


save_path = '/home/chesh/A_testchamber/working data'

# Create list of all swc morphology files in folder swc_path
swc_path = '/home/chesh/A_testchamber/Retina simulation/morphology'
morphology_list = []
for afile in os.listdir(swc_path):
    if afile.endswith(".swc"):
        swc_path_file = os.path.join(swc_path, afile)
        morphology_list.append(swc_path_file)



# Generate data files aligned with swc files
voltage_path_file = '/home/chesh/A_testchamber/simulation data/sectionValues_inc1.v'
coords_path_file = '/home/chesh/A_testchamber/simulation data/sectionCoords_inc1.txt'
voltage_list = mapswc(morphology_list, voltage_path_file, coords_path_file, save_path=save_path)

# If data is already aligend, use lines 35-40 instead of lines 29-31
# Create list of all data files in folder voltage_path
#voltage_path = '/home/chesh/A_testchamber/working_data'
#voltage_list = []
#for afile in os.listdir(swc_path):
#    if afile.endswith(".swc"):
#        voltage_path_file = os.path.join(voltage_path, afile[:-4]+'_aligned_data.v')
#        voltage_list.append(voltage_path_file)



# Convert swc(s) into vtk(s)
vtk_name='NeuroModel.vtk'
datatitle='vdata'
maxframes = 500
spherediv=8
cyldiv=5

numvtks, minV, maxV = swc2vtk(morphology_list, voltage_list, save_path=save_path, vtk_name=vtk_name, datatitle=datatitle, spherediv=spherediv, cyldiv=cyldiv, maxframes=maxframes)



# Convert vtk(s) and time data into animated gif(s)
vtk_path = os.path.join(save_path, vtk_name)
time_path_file = '/home/chesh/A_testchamber/simulation data/sectionTimes_inc1.txt'
#selectframes = '/home/chesh/A_testchamber/working_data/selectFrames_inc1.txt'
#selectframes = list(range(2000,2299+1))
selectframes=[]
colorBounds = (minV, maxV)
fps=50
startelevation=50
stepelevation=-0.15
startazimuth=130
stepazimuth=0.3

vtk2gif(vtk_path, selectframes=selectframes, datatitle=datatitle, save_path=save_path, numvtks=numvtks, colorBounds=colorBounds, maxframes=maxframes, fps=fps, time_path_file=time_path_file, startelevation=startelevation, stepelevation=stepelevation, startazimuth=startazimuth, stepazimuth=stepazimuth)



t1 = time()
print("\nTime Elapsed (s): ")
print(t1-t0)
