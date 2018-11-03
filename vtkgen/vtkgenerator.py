# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 12:37:20 2016

@author: nebula

Modified on Wed May  15 10:15:10 2018

@modder: Andres Morales
"""

import os
import math
import numpy as np
from tqdm import tqdm
from tqdm import trange

from swc2gif.vtkgen.genprimitives import GenPrimitives
from swc2gif.vtkgen.swc import Swc


class VtkGenerator():
    header_base = '''\
# vtk DataFile Version 3.0
vtkgen
ASCII
DATASET UNSTRUCTURED_GRID
'''
    volume_header_base = '''\
# vtk DataFile Version 1.0
vtkgen VOLUME
ASCII
DATASET STRUCTURED_POINTS
'''

    def __init__(self):
        self.header = self.header_base
        self.volume_header = self.volume_header_base
        
        self.swc_list = []
        self.point_list = []
        self.cell_list = []
        self.datafile_list = []
        self.timeddatafile_list = []
        self.annotation_point_list = []
        self.annotation_cell_list = []

        self.point_text = ''
        self.cell_text = ''

        self.converted = False
        self.draw_mode = 3
        self.ncell_per_compartment = []
        self.ncell_per_compartment_per_swc = []

    def set_draw_mode(self, draw_mode):
        self.draw_mode = draw_mode
        self.converted = False
        self.point_list = []
        self.cell_list = []

    def add_cylinder(self, pos_x=0, pos_y=0, pos_z=0, radius=1.0, height=1.0, rot_y=0, rot_z=0,
                     data=0.0, radius_ratio=1.0):

        local_cell_list = []
        local_point_list = []
        if self.draw_mode == 0:
            local_cell_list, local_point_list = GenPrimitives.cylinder()
            self.ncell_per_compartment.append(1)
        elif self.draw_mode == 1:
            local_cell_list, local_point_list = GenPrimitives.cylinder(top_face_diam=radius_ratio)
            self.ncell_per_compartment.append(1)
        elif self.draw_mode == 2:
            local_cell_list, local_point_list = GenPrimitives.cylinder_3cell(top_face_diam=radius_ratio)
            self.ncell_per_compartment.append(len(local_cell_list))
        elif self.draw_mode == 3:
            local_cell_list, local_point_list = GenPrimitives.hemisphere_cylinder(div=7, top_face_diam=radius_ratio,
                                                                                          height=height, radius=radius)
            self.ncell_per_compartment.append(len(local_cell_list))
            height = 1.0
            radius = 1.0
        elif self.draw_mode == 4:
            local_cell_list, local_point_list = GenPrimitives.line()
            self.ncell_per_compartment.append(1)

        point_start = len(self.point_list)
        for cell in local_cell_list:
            cell['points'] = [(i + point_start) for i in cell['points']]
            cell['data'] = data

        # scale
        local_point_list = np.array([ v*[height, radius, radius] for v in local_point_list])
        # rot z
        local_point_list = np.array([ [v[0]*np.cos(rot_z)-v[1]*np.sin(rot_z), v[0]*np.sin(rot_z)+v[1]*np.cos(rot_z), v[2]] for v in local_point_list])
        # rot y
        local_point_list = np.array([ [v[0]*np.cos(rot_y)+v[2]*np.sin(rot_y), v[1], -v[0]*np.sin(rot_y)+v[2]*np.cos(rot_y)] for v in local_point_list])
        # move
        local_point_list = np.array([ v+[pos_x, pos_y, pos_z] for v in local_point_list])

        self.point_list.extend(local_point_list)
        self.cell_list.extend(local_cell_list)

    def add_cylinder_p2p(self, pos1=(0, 0, 0), pos2=(2, 0, 0), size=1.0, data=0, radius_ratio=1.0):
        pos1 = np.array(pos1)
        pos2 = np.array(pos2)
        local_pos = pos2 - pos1

        rot_y = -np.arctan2(local_pos[2], local_pos[0])
        rot_z = np.arctan2(local_pos[1], np.sqrt(local_pos[0]**2 + local_pos[2]**2))
        length = np.sqrt(local_pos[0]**2 + local_pos[1]**2 + local_pos[2]**2)

        self.add_cylinder(pos1[0], pos1[1], pos1[2], size, length, rot_y, rot_z, data, radius_ratio=radius_ratio)

        
        
    def add_sphere(self, pos=(0,0,0), radius=1.0, data=0.0):
        pos = np.array(pos)
        
        local_cell_list = []
        local_point_list = []
            
        if self.draw_mode == 0:
            local_cell_list, local_point_list = GenPrimitives.sphere()
            self.ncell_per_compartment.append(1)
        elif self.draw_mode == 1:
            local_cell_list, local_point_list = GenPrimitives.sphere()
            self.ncell_per_compartment.append(1)
        elif self.draw_mode == 2:
            local_cell_list, local_point_list = GenPrimitives.sphere(pos, size=radius, data=data)
            self.ncell_per_compartment.append(len(local_cell_list))
        elif self.draw_mode == 3:
            local_cell_list, local_point_list = GenPrimitives.sphere(pos, size=radius, data=data)
            self.ncell_per_compartment.append(len(local_cell_list))
        elif self.draw_mode == 4:
            local_cell_list, local_point_list = GenPrimitives.sphere()
            self.ncell_per_compartment.append(1)

        point_start = len(self.point_list)
        for cell in local_cell_list:
            cell['points'] = [(i + point_start) for i in cell['points']]
            cell['data'] = data

        self.point_list.extend(local_point_list)
        self.cell_list.extend(local_cell_list)
        
        
        
    def convert_swc(self, diam_ratio=1.0, normalize_diam=False):
        self.converted = True
        self.ncell_per_compartment_per_swc = []
        
        for swc_data in tqdm(self.swc_list, desc='Convert SWC Data Into Primitive Shapes'):
            data_size = len(swc_data.data)
            
            for record in tqdm(swc_data.data.values(), desc='Converting: ' + swc_data.filename):
                if normalize_diam:
                    drawing_radius = math.sqrt(record['radius']) * diam_ratio
                else:
                    drawing_radius = record['radius'] * diam_ratio
                    
                if record['parent'] <= 0 or record['type'] == 1:
                    self.add_sphere(record['pos'], drawing_radius, float(record['id']) / data_size)
                else:
                    parent_record = swc_data.data[record['parent']]
                    if normalize_diam:
                        drawing_parent_radius = math.sqrt(parent_record['radius']) * diam_ratio
                    else:
                        drawing_parent_radius = parent_record['radius'] * diam_ratio
                    
                    self.add_cylinder_p2p(record['pos'], parent_record['pos'], drawing_radius,
                                          float(record['id']) / data_size,
                                          radius_ratio=(drawing_parent_radius/drawing_radius))
            self.ncell_per_compartment_per_swc.append(self.ncell_per_compartment)
            self.ncell_per_compartment = []
        self.point_text = self._point2text(self.point_list)
        self.cell_text = self._cell2text(self.cell_list)

    def add_swc(self, swc_filename,
                shift_x=0.0, shift_y=0.0, shift_z=0.0, inv_x=False, inv_y=False, inv_z=False):
        """add swc file for generating vtk file

        :param swc_filename: swc filename
        :param shift_x:
        :param shift_y:
        :param shift_z:
        :param inv_x:
        :param inv_y:
        :param inv_z:
        :return:
        """
        self.converted = False
        self.swc_list.append(Swc(swc_filename))
        self.swc_list[-1].invert(inv_x, inv_y, inv_z)
        self.swc_list[-1].shift(shift_x, shift_y, shift_z)

    @staticmethod
    def _point2text(point_list):
        text = 'POINTS %d float\n' % (len(point_list))
        for point in tqdm(point_list, desc='Generating Point Text'):
            text += '%f %f %f\n' % (point[0], point[1], point[2])
            
        return text

    @staticmethod
    def _cell2text(cell_list, title='data'):
        num_data = sum([len(cell['points'])+1 for cell in cell_list])
        
        text = '\nCELLS %d %d\n' % (len(cell_list), num_data)
        for cell in tqdm(cell_list, desc='Generating Cell Text'):
            text += str(len(cell['points']))
            for x in cell['points']:
                text += ' '+str(x)

            text += '\n'

        text += '\nCELL_TYPES %d\n' % (len(cell_list))
        for cell in tqdm(cell_list, desc='Generating Cell Type Text'):
            text += str(cell['type'])+'\n'

        text += '\nCELL_DATA %d\n' % (len(cell_list))
        text += 'SCALARS ' + title + ' float 1\n'
        text += 'LOOKUP_TABLE default\n'
        for cell in tqdm(cell_list, desc='Generating Cell Data Text'):
            text += str(cell['data'])+'\n'
        
        return text

    @staticmethod
    def _fixedval2text(cell_list, title='fixedval', fixedval=0.0):
        text = ''
        text += 'SCALARS '+title+' float 1\n'
        text += 'LOOKUP_TABLE default\n'
        for cell in cell_list:
            text += str(fixedval)+'\n'

        return text

    @staticmethod
    def _movingval2text(cell_list, title='movingval', num=100):
        text = ''
        text += 'SCALARS '+title+' float ' + str(num) + '\n'
        text += 'LOOKUP_TABLE default\n'

        for i in range(len(cell_list)):
            val = i
            for j in range(num):
                text += str(val) + ' '
                val = (val + 1) % 512

            text += '\n'

        return text

    def _file2text(self, datafile_list, title):
        text = ''
        text += 'SCALARS ' + title + ' float 1\n'
        text += 'LOOKUP_TABLE default\n'
        
        if not (len(datafile_list) == len(self.swc_list)):
            print('Warning: there is mismatch of data file and swc file (datafile=%d, swcfile=%d)'
                  % (len(datafile_list), len(self.swc_list)))
        for datafile_list_index, filename in enumerate(datafile_list):
            data_num = 0
            with open(filename, 'r') as f:
                read_data = f.readlines()

            for i in range(len(read_data)):
                if read_data[i][0] != '#':
                    for j in range(self.ncell_per_compartment[data_num]):
                        text += read_data[i].rstrip() + '\n'
                    data_num += 1

            if not (data_num == len(self.swc_list[datafile_list_index].data)):
                print('Warning: there is mismatch of data file lines and swc file lines (index=%d, datafile=%d, swcfile=%d)'
                      % (datafile_list_index, data_num, len(self.swc_list[datafile_list_index].data)))
        return text
        
    def _timedfile2text(self, filename, datafile_list, title, delimiter, maxframes):
        
        # Check if number of data files is 1 or equal to the number of swc's
        if len(datafile_list) != 1 and len(datafile_list) != len(self.swc_list):
            print('Warning: there is mismatch of data file and swc file (datafile=%d, swcfile=%d)'
                  % (len(datafile_list), len(self.swc_list)))
        
#        # TO DO: Error Checking!!
#        # Check if data files have same number of data points as swc files have segments   
#        # Determine number of time steps # TO DO: repeat for all files and check consistency
        with open(datafile_list[0], 'r') as f:
            read_data = f.readlines()
            timesteps = 0
            for i in range(len(read_data)):
                if read_data[i][0] != ' ':
                    timesteps += 1
            if timesteps == 0:
                print('Warning: there is no data in file (timeddatafile=%d)'
                      % (0))
        
        timesteps_digits = len(str(timesteps))
        
            
        # Calculate the number of splits
        from math import ceil
        if maxframes > timesteps:
            maxframes = timesteps
        num_vtks = int(ceil(timesteps/float(maxframes)))
        
        # Generate names of split VTK files
        vtk_path_list = []
        if num_vtks > 1:
            for splitCount in range(num_vtks):
                vtk_path_list.append(filename[:-4]+'_'+str(splitCount)+filename[-4:])
        else:
            vtk_path_list.append(filename)
        
        # Initialize variables for finding the boundary values for future colorbars
        tempMin = float('Inf')
        tempMax = -float('Inf')
        
        # Loop through each data file, generating a line offset lookup table for hopefully faster searching
        line_offset_file_list = []
        for datafilename in tqdm(datafile_list, desc='Generating Data Lookup-Table'):
            line_offset = []
            offset = 0
            with open(datafilename, 'r') as f:
                for line in f:
                    line_offset.append(offset)
                    offset += len(line) #for python 2.7 need +1 because len() will index of '\n' not the start of the new line
            line_offset_file_list.append(list(line_offset))
        
        # Loop through each VTK path, generating and writing data to the VTK
        for splitCount, vtk_path in enumerate(tqdm(vtk_path_list, desc='Writing VTK File(s)')):
            with open(vtk_path, 'w') as vf:
                
                # Write header and morphology to this vtk file
                vf.write(self.header)
                vf.write(self.point_text)
                vf.write(self.cell_text)
                
                # Calculate time bounds for this vtk file
                startstep = maxframes*splitCount
                if splitCount+1 == num_vtks:
                    stopstep = timesteps
                else:
                    stopstep = maxframes*(splitCount+1)
                
                # Iterate through each time step, writing to this vtk file
                tempHead, vtk_file = os.path.split(vtk_path)
                #for step in trange(startstep,stopstep, desc='Current [Min, Max]: ['+str(tempMin)+', '+str(tempMax)+']; Writing '+vtk_file):
                for step in trange(startstep,stopstep, desc='Writing '+vtk_file):
                    
                    # Add leading 0's for string formatting of time order for filedata name
                    stepstring = str(step+1)
                    if len(stepstring) < timesteps_digits:
                        stepstring = '0'*(timesteps_digits-len(stepstring))+stepstring
                    # Generate header
                    vf.write('SCALARS ' + title + stepstring + ' float 1\n')
                    vf.write('LOOKUP_TABLE default\n')
                    #text_list.append('SCALARS ' + title + stepstring + ' float 1')
                    #text_list.append('LOOKUP_TABLE default')
                    
                    # Iterate through each data file, reorienting and writing into a single column
                    for datafile_list_index, datafilename in enumerate(datafile_list):
                        line = ''
                        line_location = line_offset_file_list[datafile_list_index][step]
                        with open(datafilename, 'r') as f:
                            f.seek(line_location)
                            line = f.readline()
                            
                        compartments_data = line.rstrip().split(delimiter)
                        for compartment_idx, compartment_val in enumerate(compartments_data):
                            # Also check for min and max values
                            try:
                                tempMin = min(tempMin, float(compartment_val))
                                tempMax = max(tempMax, float(compartment_val))
                            except:
                                print('Data File: '+vtk_file)
                                print('Step: '+str(step))
                                print('Line offset: '+str(line_location))
                                print('Line: '+line)
                                print('Compartment Idx: '+str(compartment_idx))
                                print('Compartment Value: '+compartment_val+':')
                            for j in range(self.ncell_per_compartment_per_swc[datafile_list_index][compartment_idx]):
                                vf.write(compartment_val+'\n')
                                #text_list.append(compartment_val)
        
        #print('splitCount:')
        #print(splitCount)
        #print('Maxframes:')
        #print(maxframes)
        #print('Step:')
        #print(step)
        return [num_vtks, tempMin, tempMax]

    def _radius2text(self):
        text = ''
        text += 'SCALARS radius float 1\n'
        text += 'LOOKUP_TABLE default\n'

        for swc in self.swc_list:
            for j, record in swc.data.items():
                if j > 1:
                    for k in range(self.ncell_per_compartment[j]):
                        text += str(record['radius']) + '\n'

        return text

    def _type2text(self):
        text = ''
        text += 'SCALARS type float 1\n'
        text += 'LOOKUP_TABLE default\n'

        for swc in self.swc_list:
            for j, record in swc.data.items():
                if j > 1:
                    for k in range(self.ncell_per_compartment[j]):
                        text += str(record['type']) + '\n'

        return text

    def _coloringbyswc(self):
        text = ''
        text += 'SCALARS coloring float 1\n'
        text += 'LOOKUP_TABLE default\n'

        for i, swc in enumerate(self.swc_list):
            val = i * (1.0 / len(self.swc_list))
            for j in range(len(swc.data) - 1):
                for k in range(self.ncell_per_compartment[j]):
                    text += str(val + (1.0 / len(self.swc_list) / (len(swc.data) - 1))) + '\n'

        return text

    def add_datafile(self, datafilename):
        self.datafile_list.append(datafilename)

    def clear_datafile(self):
        self.datafile_list = []
    
    def add_timeddatafile(self, datafilename):
        self.timeddatafile_list.append(datafilename)

    def clear_timeddatafile(self):
        self.timeddatafile_list = []

    def add_mark(self, pos=(0, 0, 0), size=1.0, data=0.0):
        local_cell_list, local_point_list = \
            GenPrimitives.sphere(pos, size=size, data=data, point_start=len(self.annotation_point_list))
        self.annotation_cell_list.extend(local_cell_list)
        self.annotation_point_list.extend(local_point_list)

    def add_swc_mark(self, swc_index, compartment_index, size=1.0, data=0.0):
        if swc_index < len(self.swc_list):
            # print self.swc_list[swc_index].data
            if compartment_index in self.swc_list[swc_index].data:
                self.add_mark(self.swc_list[swc_index].data[compartment_index]['pos'], size=size, data=data)
                return self.swc_list[swc_index].data[compartment_index]['pos']
            else:
                print('Warning: Invalid compartment index (swc_id=%d, compartment_id=%d)' % (swc_index, compartment_index))
        else:
            print('Warning: Invalid swc index (swc_id=%d)' % swc_index)

    def add_swc_connection(self, swc_index1, swc_compartment1, swc_index2, swc_compartment2, size=1.0, data=1.0):
        pos1 = self.add_swc_mark(swc_index1, swc_compartment1, size, data)
        pos2 = self.add_swc_mark(swc_index2, swc_compartment2, size, data)
        local_cell_list, local_point_list = \
            GenPrimitives.line(pos1, pos2, data, point_start=len(self.annotation_point_list))

        self.annotation_cell_list.extend(local_cell_list)
        self.annotation_point_list.extend(local_point_list)

    def write_annotation_vtk(self, filename):
        with open(filename, 'w') as wfile:
            wfile.write(self.header)
            wfile.write(self._point2text(self.annotation_point_list))
            wfile.write(self._cell2text(self.annotation_cell_list, title='annotation_data'))

    def write_swc(self, filename, swc_index=0, comment='vtkgen'):
        swc = self.swc_list[swc_index]
        with open(filename, 'w') as swcfile:
            swcfile.write(swc.header)
            swcfile.write('# ' + comment + '\n')
            for j, record in swc.data.items():
                swcfile.write('%d %d %f %f %f %f %d\n' % (record['id'], record['type'],
                                                        record['pos'][0], record['pos'][1], record['pos'][2],
                                                        record['radius'], record['parent']))

    def write_vtk(self, filename, fixedval=None, datatitle='filedata', movingval=False, coloring=False,
                  diam_ratio=1.0, normalize_diam=False, radius_data=False, type_data=False, delimiter=' ',
                  maxframes=float('Inf')):
        """generate and write vtk to file

        :param filename: Output VTK filename
        :param fixedval: add fixed value to CELL_DATA in VTK file
        :param datatitle: change title of CELL_DATA from appended data file
        :param movingval: add moving value to CELL_DATA in VTK file
        :param coloring: add coloring data to CELL_DATA in VTK file
        :param diam_ratio: multiply diam_ratio to all diameter of SWC compartments
        :param normalize_diam: sqrt(diam) to normalize diameter of SWC compartments
        :param radius_data: add radius information of SWC compartments to CELL_DATA in VTK file
        :param type_data: Output type information of SWC compartments to CELL_DATA in VTK file
        :type filename: text
        :type fixedval: float
        :type datatitle: text
        :type movingval: bool
        :type coloring: bool
        :type diam_ratio: float
        :type normalize_diam: bool
        :type radius_data: bool
        :type type_data: bool
        :return:
        """
        if not self.converted:
            self.convert_swc(diam_ratio=diam_ratio, normalize_diam=normalize_diam)
        
        if len(self.timeddatafile_list) > 0:
            num_vtks, minV, maxV = self._timedfile2text(filename, self.timeddatafile_list, datatitle, delimiter, maxframes)
            
            return [num_vtks, minV, maxV]
        else:
            with open(filename, 'w') as file:
                file.write(self.header)
                file.write(self.point_text)
                file.write(self.cell_text)

                if fixedval is not None:
                    file.write(self._fixedval2text(self.cell_list, fixedval=fixedval))

                if movingval:
                    file.write(self._movingval2text(self.cell_list))

                if len(self.datafile_list) > 0:
                    file.write(self._file2text(self.datafile_list, datatitle))

                if coloring:
                    file.write(self._coloringbyswc())

                if radius_data:
                    file.write(self._radius2text())

                if type_data:
                    file.write(self._type2text())
            return [1, minV, maxV]

    @staticmethod
    def _swc2volume(swc, world, origin=(0.0, 0.0, 0.0), ratio=(1.0, 1.0, 1.0), point_weight=0.2):
        for k, record in tqdm(swc.data.items(), desc=swc.filename):
            pos = (int(round((record['pos'][0] - origin[0]) / ratio[0], 0)),
                   int(round((record['pos'][1] - origin[1]) / ratio[1], 0)),
                   int(round((record['pos'][2] - origin[2]) / ratio[2], 0)))
            if pos[2] < 0 or pos[1] < 0 or pos[0] < 0\
                    or pos[2] >= len(world) or pos[1] >= len(world[0]) or pos[0] >= len(world[0][0]):
                print('Out of range: (%f, %f, %f)' % (record['pos'][0], record['pos'][1], record['pos'][2]))
            else:
                world[pos[2]][pos[1]][pos[0]] += point_weight

        return world

    def write_volume_vtk(self, filename, origin=(0.0, 0.0, 0.0), ratio=(1.0, 1.0, 1.0), div=(256, 256, 64)):
        vtkdata = ''
        vtkdata += self.volume_header
        vtkdata += 'DIMENSIONS %d %d %d\n' % (div[0], div[1], div[2])
        vtkdata += 'ORIGIN %f %f %f\n' % (origin[0], origin[1], origin[2])
        vtkdata += 'ASPECT_RATIO %f %f %f\n' % (ratio[0], ratio[1], ratio[2])
        vtkdata += ('POINT_DATA %d\n' % (div[0] * div[1] * div[2]))
        vtkdata += 'SCALARS volume float\n'
        vtkdata += 'LOOKUP_TABLE default\n'

        world = np.zeros((div[2], div[1], div[0]))

        for swc in self.swc_list:
            world = self._swc2volume(swc, world, origin, ratio)

        with open (filename, 'w') as file:
            file.write(vtkdata)
            for i in range(div[2]):
                for j in range(div[1]):
                    for k in range(div[0]):
                        file.write('%f\n' % (world[i][j][k]))

    def show_state(self):
        print(self.point_list)
        print(self.cell_list)
        print(self.swc_list)
        print(self.datafile_list)


if __name__ == '__main__':

    def test_swc_movie(stoptime=100):
        filename_base = 'swc_cuboid%d.vtk'
        vtkgen = VtkGenerator()
        
        add_swc(os.path.join('swc', 'Swc_BN_1056.swc'))
    
        for t in range(stoptime):
            print('t = %d' % t)
            for i in range(len(cell_list)):
                cell_list[i]['data'] += 0.01
                if cell_list[i]['data'] > 1.0:
                    cell_list[i]['data'] = 0.0
                
            write_vtk(filename_base % t)
    
    
    def test_cylinder(num):
        filename = 'cylinder.vtk'
        vtkgen = VtkGenerator()

        for i in range(num):
            add_cylinder(i*5, 0, 0, i*0.1+1, i*2+1, 0, 0, 0.2*i)
        
        write_vtk(filename)
        
    test_cylinder(1)
    # test_swc_movie(100)
    # test_swc_line()
