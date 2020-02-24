#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 Zhou Dejun
# The funcitons aim at following purposes.
#1.Read DICOM-CT files, DICOM-structure file, DICOM-dose-storage, 
#  DICOM-dose-plan file.
#2.Input [yaw, pitch, yaw] rotation values 
#3.Based on the values, rotate the 3D volume around a selected point,
#  typically the isodose center.
#4.Export 2D files
#5.The export files can be imported into treatment planning system
"""
Last modified on 2020-2-24

@author: Zhou Dejun
"""

import os
import numpy as np
import pydicom
from pyquaternion import Quaternion
from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

class ArrayVolume(object):
     def __init__(self, Origin = None, Pivot= None,\
                  Vector_X= None, Vector_Y= None, Vector_Z = None,\
                  Col_Spacing= None, Row_Spacing= None, Thickness = None,\
                  Array_Volume = None, Data_Type = None):
        self.Origin = Origin
        self.Pivot = Pivot
        self.Vector_X = Vector_X
        self.Vector_Y = Vector_Y
        self.Vector_Z = Vector_Z
        self.Col_Spacing = Col_Spacing
        self.Row_Spacing = Row_Spacing
        self.Thickness = Thickness
        self.Array_Volume= Array_Volume
        self.Data_Type = Data_Type
        
def build_image_volume(folder_path):
    image_volume = ArrayVolume()
    file_list = finder_ct_image_file(folder_path)
    init_image_ds= pydicom.dcmread(file_list[0])
    col_spacing = init_image_ds.PixelSpacing[1]
    row_spacing = init_image_ds.PixelSpacing[0]
    slice_thick = init_image_ds.SliceThickness
    start_position = pydicom.dcmread(file_list[0]).ImagePositionPatient
    origin = np.array(start_position)-\
              np.array([col_spacing/2, row_spacing/2, slice_thick/2])
    pivot = pydicom.dcmread(finder_rt_plan_file(folder_path))\
            [0x300a,0xb0][0][0x300a,0x111][0][0x300a,0x12c].value
    vector_x = [1,0,0]
    vector_y = [0,1,0]
    vector_z = [0,0,1]
    data_type = 'uint'+str(init_image_ds.BitsStored)
    image_volume_array = []
    for frame in range(0, len(file_list)):
        image_volume_array.append\
            (pydicom.dcmread(file_list[frame]).pixel_array)
    image_volume_array = np.array(image_volume_array)
    
    image_volume = ArrayVolume(origin, pivot, vector_x, vector_y, vector_z,\
                        col_spacing, row_spacing, slice_thick,
                        image_volume_array, data_type)  

    return image_volume  

def build_dose_volume(folder_path):
    dose_volume = ArrayVolume()
    ds= pydicom.dcmread(finder_rt_dose_file(folder_path))
    col_spacing = ds.PixelSpacing[1]
    row_spacing = ds.PixelSpacing[0]
    slice_thick = ds.GridFrameOffsetVector[1]-ds.GridFrameOffsetVector[0]
    start_position = ds.ImagePositionPatient
    origin = np.array(start_position)-\
              np.array([col_spacing/2, row_spacing/2, slice_thick/2])
    pivot = pydicom.dcmread(finder_rt_plan_file(folder_path))\
            [0x300a,0xb0][0][0x300a,0x111][0][0x300a,0x12c].value
    vector_x = [1,0,0]
    vector_y = [0,1,0]
    vector_z = [0,0,1]
    dose_volume_array = np.array(ds.pixel_array)
    data_type = 'uint'+str(ds.BitsStored)
    dose_volume = ArrayVolume(origin, pivot, vector_x, vector_y, vector_z,\
                        col_spacing, row_spacing, slice_thick,
                        dose_volume_array, data_type)  

    return dose_volume

def rotate_array_volume(array_volume, yaw, pitch, roll):
    rotated_array_volume = ArrayVolume()
    yaw = np.deg2rad(360-yaw)
    pitch = np.deg2rad(pitch)
    roll = np.deg2rad(360-roll)
    pitch_rotate = Quaternion(axis=[1, 0, 0], angle=pitch)
    roll_rotate = Quaternion(axis=[0, 1, 0], angle=roll)
    yaw_rotate = Quaternion(axis=[0, 0, 1], angle=yaw)
    rotation = pitch_rotate * roll_rotate * yaw_rotate  
    vector_origin = np.array(array_volume.Origin)-np.array(array_volume.Pivot)
    vector_origin_rot = rotation.rotate(vector_origin)
    rotated_array_volume.Origin= (np.array(vector_origin_rot)+\
                  np.array(array_volume.Pivot)).tolist()
    rotated_array_volume.Vector_X = rotation.rotate(array_volume.Vector_X)
    rotated_array_volume.Vector_Y = rotation.rotate(array_volume.Vector_Y)
    rotated_array_volume.Vector_Z = rotation.rotate(array_volume.Vector_Z)
    rotated_array_volume.Pivot = array_volume.Pivot
    rotated_array_volume.Col_Spacing = array_volume.Col_Spacing
    rotated_array_volume.Row_Spacing = array_volume.Row_Spacing
    rotated_array_volume.Thickness = array_volume.Thickness
    rotated_array_volume.Array_Volume = array_volume.Array_Volume
    rotated_array_volume.Data_Type = array_volume.Data_Type
      
    return rotated_array_volume

def map_rotate_array_volume(array_volume):       
    array_map_frames = array_volume.Array_Volume.shape[0]
    array_map_rows = array_volume.Array_Volume.shape[1]
    array_map_cols = array_volume.Array_Volume.shape[2]
    col_spacing = array_volume.Col_Spacing
    row_spacing = array_volume.Row_Spacing
    thickness = array_volume.Thickness
    vector_x = array_volume.Vector_X
    vector_y = array_volume.Vector_Y
    vector_z = array_volume.Vector_Z
    init_position = np.array(array_volume.Origin)+\
                    (col_spacing/2)*np.array(vector_x)+\
                    (row_spacing/2)*np.array(vector_y)+\
                    (thickness/2)*np.array(vector_z)    
    array_map = np.zeros([array_map_frames, array_map_rows, array_map_cols, 3])
    for frames in range(0, array_map_frames):
        print('Finished mapping '+str(frames+1)+' frames. Total '+\
              str(array_map_frames)+' .')
        for rows in range(0, array_map_rows):
            for cols in range(0, array_map_cols): 
                array_map[frames][rows][cols] = init_position+\
                    (col_spacing*cols)*np.array(vector_x)+\
                    (row_spacing*rows)*np.array(vector_y)+\
                    (thickness*frames)*np.array(vector_z)
    
    return array_map
  
def resample_rotate_array_volume_zplane(array_volume, yaw, pitch, roll):
    rotated_volume = rotate_array_volume(array_volume, yaw, pitch, roll)
    map_rotated_volume = map_rotate_array_volume(rotated_volume)
    vertex = volume_vertex(rotated_volume)
    x_list = []
    y_list = []
    z_list = []
    for i in range(0, len(vertex)):
        x_list.append(vertex[i][0])
        y_list.append(vertex[i][1])
    for i in range(0, array_volume.Array_Volume.shape[0]):
        z_list.append(array_volume.Origin[2]+array_volume.Thickness*i)
    x_min = np.min(x_list)
    y_min = np.min(y_list)
    z_min = np.min(z_list)
    z_max = np.max(z_list)+array_volume.Thickness
    slice_cols = int((np.max(x_list)-x_min)/array_volume.Col_Spacing+1)
    slice_rows = int((np.max(y_list)-y_min)/array_volume.Row_Spacing+1)
    slice_frames = int((np.max(z_list)-z_min)/array_volume.Thickness+1)   
    resampled_array =np.zeros([slice_frames, slice_rows, slice_cols],\
                              dtype= array_volume.Data_Type)
    for frames in range(0, map_rotated_volume.shape[0]):
        print('Finished resample '+str(frames+1)+' frames. Total '+\
              str(map_rotated_volume.shape[0])+' .')
        for rows in range(0, map_rotated_volume.shape[1]):
            for cols in range(0, map_rotated_volume.shape[2]):
                rotated_array_point = map_rotated_volume[frames][rows][cols]
                if z_min < rotated_array_point[2] < z_max:
                    x_index = int((rotated_array_point[0]-x_min)/\
                                  array_volume.Col_Spacing)
                    y_index = int((rotated_array_point[1]-y_min)/\
                                  array_volume.Row_Spacing)
                    z_index = int((rotated_array_point[2]-z_min)/\
                                  array_volume.Thickness)
                    resampled_array[z_index][y_index][x_index]=\
                        rotated_volume.Array_Volume[frames][rows][cols]
    resampled_array_volume = ArrayVolume()
    resampled_array_volume.Array_Volume = resampled_array
    resampled_array_volume.Origin = [x_min,y_min,z_min]
    resampled_array_volume.Pivot = array_volume.Pivot
    resampled_array_volume.Vector_X = [1,0,0]
    resampled_array_volume.Vector_Y = [0,1,0]
    resampled_array_volume.Vector_Z = [0,0,1]
    resampled_array_volume.Col_Spacing= array_volume.Col_Spacing
    resampled_array_volume.Row_Spacing= array_volume.Row_Spacing
    resampled_array_volume.Thickness= array_volume.Thickness
    
    return resampled_array_volume

def export_rotated_image_files(folder_path, yaw, pitch, roll, \
                               removedot = False):
    origin_file_list = finder_ct_image_file(folder_path)
    rotated_file_path = folder_path+\
                        '/'+'Y'+str(yaw)+'P'+str(pitch)+'R'+str(roll)
    if not os.path.exists(rotated_file_path):
        os.makedirs(rotated_file_path)        
    image_volume = build_image_volume(folder_path)
    rotated_image_volume =\
        resample_rotate_array_volume_zplane(image_volume, \
                                            yaw, pitch, roll)
    for files_index in range(0, len(origin_file_list)):
        ds = pydicom.dcmread(origin_file_list[files_index])
        ds.Rows = rotated_image_volume.Array_Volume.shape[1]
        ds.Columns = rotated_image_volume.Array_Volume.shape[2]
        ds.ImagePositionPatient = (np.array(rotated_image_volume.Origin)+\
                                np.array([ds.PixelSpacing[1]/2,\
                                ds.PixelSpacing[0]/2,\
                                ds.SliceThickness*(files_index+1/2)])).tolist()
        rotated_image_array = rotated_image_volume.Array_Volume[files_index]
        if removedot==True:
            rotated_image_array = remove_2d_black_dot(rotated_image_array)
        ds.PixelData = np.array(rotated_image_array).tobytes()   
        rotated_file_name = rotated_file_path+\
                      origin_file_list[files_index][len(folder_path):-4]+\
                      'Y'+str(yaw)+'P'+str(pitch)+'R'+str(roll)+'.dcm'
        ds.save_as(rotated_file_name)
        print('Export '+str(int(files_index+1))+' files of '+\
              str(len(origin_file_list))+\
              ' files in total.')
    return print('Rotated '+str(len(origin_file_list))+' ct images.')

def export_rotated_dose_files(folder_path, yaw, pitch, roll,\
                              removedot = False):
    origin_file_path = finder_rt_dose_file(folder_path)
    rotated_file_path = folder_path+\
                        '/'+'Y'+str(yaw)+'P'+str(pitch)+'R'+str(roll)    
    if not os.path.exists(rotated_file_path):
        os.makedirs(rotated_file_path)        
    dose_volume = build_dose_volume(folder_path)
    rotated_dose_volume =\
        resample_rotate_array_volume_zplane(dose_volume, \
                                            yaw, pitch, roll)
    ds = pydicom.dcmread(origin_file_path)
    frame_thickness = ds.GridFrameOffsetVector[1]-ds.GridFrameOffsetVector[0]
    ds.Rows = rotated_dose_volume.Array_Volume.shape[1]
    ds.Columns = rotated_dose_volume.Array_Volume.shape[2]
    ds.ImagePositionPatient = (np.array(rotated_dose_volume.Origin)+\
                                np.array([ds.PixelSpacing[1]/2,\
                                ds.PixelSpacing[0]/2,\
                                frame_thickness/2])).tolist()    
    if removedot == True:
        for frame in range(0, rotated_dose_volume.Array_Volume.shape[0]):
            rotated_dose_volume.Array_Volume[frame] = \
                remove_2d_black_dot(rotated_dose_volume.Array_Volume[frame])
    ds.PixelData = np.array(rotated_dose_volume.Array_Volume).tobytes()   
    rotated_file_name = rotated_file_path+\
                      origin_file_path[len(folder_path):-4]+\
                      'Y'+str(yaw)+'P'+str(pitch)+'R'+str(roll)+'.dcm'
    ds.save_as(rotated_file_name)
    dose_plan_file = finder_rt_plan_file(folder_path)
    ds = pydicom.dcmread(dose_plan_file)
    plan_file_name = rotated_file_path+\
                      dose_plan_file[len(folder_path):-4]+\
                      'Y'+str(yaw)+'P'+str(pitch)+'R'+str(roll)+'.dcm'
    ds.save_as(plan_file_name)
    return print('Rotated rt dose file and copied dose plan file.')

def volume_vertex(array_volume):
    origin = array_volume.Origin
    x_length = array_volume.Col_Spacing*array_volume.Array_Volume.shape[2]
    y_length = array_volume.Row_Spacing*array_volume.Array_Volume.shape[1]
    z_length = array_volume.Thickness*array_volume.Array_Volume.shape[0]
    vector_x = array_volume.Vector_X
    vector_y = array_volume.Vector_Y
    vector_z = array_volume.Vector_Z
    x_termin = np.array(origin)+x_length*np.array(vector_x)
    y_termin = np.array(origin)+y_length*np.array(vector_y)
    z_termin = np.array(origin)+z_length*np.array(vector_z)
    xy_tip = y_termin+x_length*np.array(vector_x)
    xz_tip = z_termin+x_length*np.array(vector_x)
    yz_tip = z_termin+y_length*np.array(vector_y)
    origin_opposite = yz_tip+x_length*np.array(vector_x)
    
    return origin, x_termin, y_termin, z_termin,\
            xy_tip, xz_tip, yz_tip, origin_opposite

def finder_ct_image_file(folder_path):
    """CT Image Storage files Finder

        Locate the CT Image Storage files in the patient folder.

        Args:
            folder_path: The path to the patient folder.

        Returns:
            Return a list of the path to the CT Image Storage files.
    """
    file_list = []
    ct_image_file_list = []
    for filename in os.listdir(folder_path):
        filepath = os.path.join(folder_path, filename)
        if filepath.endswith('.dcm'):
            file_list.append(filepath)
    for item in file_list:
        if pydicom.dcmread(item).SOPClassUID.name == 'CT Image Storage':
            ct_image_file_list.append(item)
    ct_image_file_list = sorted(ct_image_file_list, 
                                key=lambda s:\
                                pydicom.dcmread(s).ImagePositionPatient[2])
    return ct_image_file_list

def finder_rt_plan_file(folder_path):
    """RT Plan Storage file Finder

        Locate the RT Structure Set Storage File in the patient folder.

        Args:
            folder_path: The path to the patient folder.

        Returns:
            If there is no dose plan file in the folder, returns a string. 
            Else return the path to the RT Plan Storage file.
    """
    file_list = []
    plan_file_path = 'No Plan File'
    for filename in os.listdir(folder_path):
        filepath = os.path.join(folder_path, filename)
        if filepath.endswith('.dcm'):
            file_list.append(filepath)
    for item in file_list:
        if pydicom.dcmread(item).SOPClassUID.name == 'RT Plan Storage':
            plan_file_path = item

    return plan_file_path

def finder_rt_structure_file(folder_path):
    """RT Structure Set Storage file Finder

        Locate the RT Structure Set Storage File in the patient folder.

        Args:
            folder_path: The path to the patient folder.

        Returns:
            If there is no structure file in the folder, returns a string. 
            Else return the path to the RT Structure Set Storage file.
    """
    file_list = []
    structure_file_path = 'No Structure File'
    for filename in os.listdir(folder_path):
        filepath = os.path.join(folder_path, filename)
        if filepath.endswith('.dcm'):
            file_list.append(filepath)
    for item in file_list:
        if pydicom.dcmread(item).SOPClassUID.name == \
            'RT Structure Set Storage':
            structure_file_path = item

    return structure_file_path

def finder_rt_dose_file(folder_path):
    """RT Dose Storage file Finder

        Locate the RT Structure Set Storage File in the patient folder.

        Args:
            folder_path: The path to the patient folder.

        Returns:
            If there is no dose storage file in the folder, returns a string. 
            Else return the path to the RT Dose Storage file.
    """
    file_list = []
    dose_file_path = 'No Dose File'
    for filename in os.listdir(folder_path):
        filepath = os.path.join(folder_path, filename)
        if filepath.endswith('.dcm'):
            file_list.append(filepath)
    for item in file_list:
        if pydicom.dcmread(item).SOPClassUID.name == 'RT Dose Storage':
            dose_file_path = item

    return dose_file_path
 
def plot_axial_plane(folder_path, file_index, \
                     showstructure = False, showisodose = False,\
                     removedot = False, colormap = cm.bone):
    file_list  = finder_ct_image_file(folder_path)
    ds = pydicom.dcmread(file_list[file_index-1])
    x_min = ds.ImagePositionPatient[0]
    x_max = ds.ImagePositionPatient[0]+(ds.Columns+1)*ds.PixelSpacing[1]    
    y_min = ds.ImagePositionPatient[1]
    y_max = ds.ImagePositionPatient[1]+(ds.Rows+1)*ds.PixelSpacing[0] 
    z_position = ds.ImagePositionPatient[2]
    

    fig = Figure()
    ax = fig.add_subplot(111)
    ax.set_title('Slice Position in SI direction ' + str(z_position)+' cm')
    secaxx = ax.secondary_xaxis('top')
    secaxy = ax.secondary_yaxis('right')
    x_major_ticks = np.linspace(-(-x_min//100)*100, (x_max//100)*100,\
                                int((x_max//100)+(-x_min//100)+1))
    y_major_ticks = np.linspace(-(-y_min//100)*100, (y_max//100)*100,\
                                int((y_max//100)+(-y_min//100)+1))
    ax.set_xticks(x_major_ticks)
    ax.set_yticks(y_major_ticks)
    ax.set_xlabel('Posterior Direction [cm]')
    ax.set_ylabel('Right Direction [cm]')    
    secaxx.set_xlabel('Anterior Direction [cm]')
    secaxy.set_ylabel('Left Direction [cm]')
    image_array = ds.pixel_array
    if removedot == True:
        image_array = remove_2d_black_dot(ds.pixel_array)
    ax.imshow(image_array,\
              extent = (x_min, x_max, y_max, y_min),cmap = cm.bone, zorder=1)
    if showisodose == True:
        dose_storage_file = pydicom.dcmread(finder_rt_dose_file(folder_path))
        dose_plan_file = pydicom.dcmread(finder_rt_plan_file(folder_path))
        dose_storage_position = dose_storage_file.ImagePositionPatient
        dose_center = dose_plan_file\
            [0x300a,0xb0][0][0x300a,0x111][0][0x300a,0x12c].value
        dose_thickness = dose_storage_file.GridFrameOffsetVector[1]-\
            dose_storage_file.GridFrameOffsetVector[0]
        dose_z_index = int((z_position-\
                            dose_storage_position[2])/\
                                dose_thickness)
        if dose_z_index <0 or dose_z_index >= dose_storage_file.NumberOfFrames:
            print('No dose stroage for this slice.')
            showisodose = False
    if showisodose == True:
        dose_x_min = dose_storage_file.ImagePositionPatient[0]-\
            dose_storage_file.PixelSpacing[1]/2
        dose_x_max = dose_x_min+\
            (dose_storage_file.Columns+1/2)*dose_storage_file.PixelSpacing[1] 
        dose_y_min = dose_storage_file.ImagePositionPatient[1]-\
            dose_storage_file.PixelSpacing[0]/2
        dose_y_max = dose_storage_file.ImagePositionPatient[1]+\
            (dose_storage_file.Rows+1/2)*dose_storage_file.PixelSpacing[0]
        dose_array = dose_storage_file.pixel_array[dose_z_index]
        if removedot == True:
            dose_array = remove_2d_black_dot(dose_array)
        dose_distribution = dose_array*dose_storage_file.DoseGridScaling
        dose_legend = legend_dose_storage(folder_path)
        ax.contour(dose_distribution, dose_legend['dose value'],\
                    colors=dose_legend['dose color'], origin='image',\
                    linewidths = 0.5, zorder=2,\
                    extent=(dose_x_min, dose_x_max, dose_y_max, dose_y_min))
        dose_list = []
        for item in range(0, len(dose_legend['dose name'])):
            dose_list.append(Line2D([0], [0],\
                                    color=dose_legend['dose color'][item],\
                                    lw=1,\
                                    label=dose_legend['dose name'][item]))
        ax.add_artist(ax.legend(handles=dose_list, \
            bbox_to_anchor=(1.25, 1), loc='upper left', borderaxespad=0.))
        if dose_center[2] == z_position:
            ax.scatter(dose_center[0], dose_center[1], c = 'r', marker = '+',\
                        zorder=4)
            dose_center_legend = \
                [Line2D([0], [0], color='r', ls='None', marker='+', \
                        label='Isodose Center')]
            ax.add_artist(ax.legend(handles=dose_center_legend, \
                          bbox_to_anchor=(1.25, 0.35), loc='upper left',\
                          borderaxespad=0.))
    if showstructure == True:
        structure_file = pydicom.dcmread(finder_rt_structure_file(folder_path))
        structure_legend = legend_structure(structure_file)
        for name in range(0, len(structure_file.ROIContourSequence)):
            print('Loading ROI num ', name)
            if hasattr(structure_file.ROIContourSequence[name],\
                        'ContourSequence') == True:
                for ctr in range(0,\
                                  len(structure_file.
                                      ROIContourSequence[name].\
                                      ContourSequence)):
                    ctr_x = []
                    ctr_y = []
                    if structure_file.ROIContourSequence[name].\
                            ContourSequence[ctr].ContourData[2] == z_position:
                        for coord in range(0,\
                                int(len(structure_file.\
                                        ROIContourSequence[name].\
                                        ContourSequence[ctr].\
                                        ContourData)/3)):                       
                            ctr_x.append(structure_file.\
                                          ROIContourSequence[name].\
                                          ContourSequence[ctr].\
                                          ContourData[0+3*coord])
                            ctr_y.append(structure_file.\
                                          ROIContourSequence[name].\
                                          ContourSequence[ctr].\
                                          ContourData[1+3*coord])
                    if len(ctr_x)!=0:
                        ax.fill(ctr_x,ctr_y, 
                            edgecolor=structure_legend['edge color'][name],\
                            facecolor=structure_legend['face color'][name],\
                            linewidth = 1,\
                            zorder=3)
        structure_list = [] 
        for item in range(0, len(structure_file.ROIContourSequence)):
            structure_list.append(Patch(\
                            facecolor=structure_legend['edge color'][item],\
                            edgecolor='black', 
                            label=structure_legend['structure name'][item]))
        ax.add_artist(ax.legend(handles=structure_list,\
                      bbox_to_anchor=(-0.25, 1), loc='upper right', \
                      borderaxespad=0.))
           
    return fig
            
def plot_sagittal_plane(folder_path, frame_index, \
                      showstructure = False, showisodose = False,\
                      removedot = False, color_map = cm.bone):
    image_volume  = build_image_volume(folder_path)
    frame_array = image_volume.Array_Volume[:,:,frame_index-1]

    x_min = image_volume.Origin[1]
    x_max = image_volume.Origin[1]+\
        (image_volume.Array_Volume.shape[1]+1)*image_volume.Row_Spacing  
    y_min = image_volume.Origin[2]
    y_max = image_volume.Origin[2]+\
        (image_volume.Array_Volume.shape[0]+1)*image_volume.Thickness
        
    z_position = round(image_volume.Origin[0]+(frame_index-1/2)*\
        image_volume.Col_Spacing,4)  
    
    fig = Figure()
    ax = fig.add_subplot(111)
    ax.set_title('Slice Position X = ' + str(z_position)+' cm')
    secaxx = ax.secondary_xaxis('top')
    secaxy = ax.secondary_yaxis('right')
    x_major_ticks = np.linspace(-(-x_min//100)*100, (x_max//100)*100,\
                                int((x_max//100)+(-x_min//100)+1))
    y_major_ticks = np.linspace(-(-y_min//100)*100, (y_max//100)*100,\
                                int((y_max//100)+(-y_min//100)+1))
    ax.set_xticks(x_major_ticks)
    ax.set_yticks(y_major_ticks)
    ax.set_xlabel('Inferior Direction [cm]')
    ax.set_ylabel('Anteriot Direction [cm]')    
    secaxx.set_xlabel('Superior Direction [cm]')
    secaxy.set_ylabel('Posterior Direction [cm]')
    if removedot == True:
        frame_array = remove_2d_black_dot(frame_array)
    ax.imshow(frame_array , cmap=color_map, \
              extent = (x_min, x_max, y_min, y_max),\
              origin = 'lower',\
              zorder=1)
    if showisodose == True:
        dose_storage_file = pydicom.dcmread(finder_rt_dose_file(folder_path))
        dose_plan_file = pydicom.dcmread(finder_rt_plan_file(folder_path))
        dose_storage_position = dose_storage_file.ImagePositionPatient
        dose_thickness = dose_storage_file.GridFrameOffsetVector[1]-\
            dose_storage_file.GridFrameOffsetVector[0]
        dose_center = dose_plan_file\
            [0x300a,0xb0][0][0x300a,0x111][0][0x300a,0x12c].value
        dose_z_index = int((z_position-dose_storage_position[0])/\
                                dose_storage_file.PixelSpacing[1])
        if dose_z_index <0 or dose_z_index >= dose_storage_file.Columns:
            print('No dose stroage for this slice.')
            showisodose = False
    if showisodose == True:
        dose_x_min = dose_storage_position[1]-\
            dose_storage_file.PixelSpacing[1]/2
        dose_x_max = dose_x_min + \
            (dose_storage_file.Rows+1)*dose_storage_file.PixelSpacing[0]
        dose_y_min = dose_storage_position[2]-dose_thickness/2
        dose_y_max = dose_y_min+\
            (dose_storage_file.NumberOfFrames+1)*dose_thickness
        dose_array = dose_storage_file.pixel_array[:,:,dose_z_index]
        if removedot == True:
            dose_array = remove_2d_black_dot(dose_array)
        dose_distribution = dose_array*dose_storage_file.DoseGridScaling
        dose_legend = legend_dose_storage(folder_path)
        ax.contour(dose_distribution, dose_legend['dose value'],\
                    colors=dose_legend['dose color'], origin='image',\
                    linewidths = 0.5, zorder=2,\
                    extent=(dose_x_min, dose_x_max, dose_y_max, dose_y_min))
        dose_list = []
        for item in range(0, len(dose_legend['dose name'])):
            dose_list.append(Line2D([0], [0],\
                                    color=dose_legend['dose color'][item],\
                                    lw=1,\
                                    label=dose_legend['dose name'][item]))
        ax.add_artist(ax.legend(handles=dose_list, \
            bbox_to_anchor=(1.25, 1), loc='upper left', borderaxespad=0.))
        if int(dose_center[0]) == int(z_position):
            ax.scatter(dose_center[1], dose_center[2], c = 'r', marker = '+',\
                        zorder=4)
            dose_center_legend = \
                [Line2D([0], [0], color='r', ls='None', marker='+', \
                        label='Isodose Center')]
            ax.add_artist(ax.legend(handles=dose_center_legend, \
                    bbox_to_anchor=(1.25, 0.35), loc='upper left', \
                    borderaxespad=0.))
    return fig

def plot_coronal_plane(folder_path, frame_index, \
                      showstructure = False, showisodose = False,
                      removedot = False, colormap = cm.bone):
    image_volume  = build_image_volume(folder_path)
    frame_array = image_volume.Array_Volume[:,frame_index-1,:]

    x_min = image_volume.Origin[0]
    x_max = image_volume.Origin[0]+\
        (image_volume.Array_Volume.shape[2]+1)*image_volume.Col_Spacing  
    y_min = image_volume.Origin[2]
    y_max = image_volume.Origin[2]+\
        (image_volume.Array_Volume.shape[0]+1)*image_volume.Thickness
        
    z_position = round(image_volume.Origin[1]+(frame_index-1/2)*\
        image_volume.Row_Spacing,4)  
    
    fig = Figure()
    ax = fig.add_subplot(111)
    ax.set_title('Slice Position Y = ' + str(z_position)+' cm')
    secaxx = ax.secondary_xaxis('top')
    secaxy = ax.secondary_yaxis('right')
    x_major_ticks = np.linspace(-(-x_min//100)*100, (x_max//100)*100,\
                                int((x_max//100)+(-x_min//100)+1))
    y_major_ticks = np.linspace(-(-y_min//100)*100, (y_max//100)*100,\
                                int((y_max//100)+(-y_min//100)+1))
    ax.set_xticks(x_major_ticks)
    ax.set_yticks(y_major_ticks)
    ax.set_xlabel('Inferior Direction [cm]')
    ax.set_ylabel('Right Direction [cm]')    
    secaxx.set_xlabel('Superior Direction [cm]')
    secaxy.set_ylabel('Left Direction [cm]')
    if removedot == True:
        frame_array = remove_2d_black_dot(frame_array)
    ax.imshow(frame_array, cmap=colormap, \
              extent = (x_min, x_max, y_min, y_max),\
              origin = 'lower',\
              zorder=1)
    if showisodose == True:
        dose_storage_file = pydicom.dcmread(finder_rt_dose_file(folder_path))
        dose_plan_file = pydicom.dcmread(finder_rt_plan_file(folder_path))
        dose_storage_position = dose_storage_file.ImagePositionPatient
        dose_thickness = dose_storage_file.GridFrameOffsetVector[1]-\
            dose_storage_file.GridFrameOffsetVector[0]
        dose_center = dose_plan_file\
            [0x300a,0xb0][0][0x300a,0x111][0][0x300a,0x12c].value
        dose_z_index = int((z_position-\
                            (dose_storage_position[1])/\
                                dose_storage_file.PixelSpacing[0]))
        if dose_z_index <0 or dose_z_index >= dose_storage_file.Rows:
            print('No dose stroage for this slice.')
            showisodose = False
    if showisodose == True:
        dose_x_min = dose_storage_position[0]-\
            dose_storage_file.PixelSpacing[1]/2
        dose_x_max = dose_x_min + \
            (dose_storage_file.Columns+1)*dose_storage_file.PixelSpacing[1]
        dose_y_min = dose_storage_position[2]-dose_thickness/2
        dose_y_max = dose_y_min+\
            (dose_storage_file.NumberOfFrames+1)*dose_thickness
        dose_array = dose_storage_file.pixel_array[:,dose_z_index,:]
        if removedot == True:
            dose_array = remove_2d_black_dot(dose_array)
        dose_distribution = dose_array*dose_storage_file.DoseGridScaling
        dose_legend = legend_dose_storage(folder_path)
        ax.contour(dose_distribution, dose_legend['dose value'],\
                    colors=dose_legend['dose color'], origin='image',\
                    linewidths = 0.5, zorder=2,\
                    extent=(dose_x_min, dose_x_max, dose_y_max, dose_y_min))
        dose_list = []
        for item in range(0, len(dose_legend['dose name'])):
            dose_list.append(Line2D([0], [0],\
                                    color=dose_legend['dose color'][item],\
                                    lw=1,\
                                    label=dose_legend['dose name'][item]))
        ax.add_artist(ax.legend(handles=dose_list, \
            bbox_to_anchor=(1.25, 1), loc='upper left', borderaxespad=0.))
        if int(dose_center[1]) == int(z_position):
            ax.scatter(dose_center[0], dose_center[2], c = 'r', marker = '+',\
                        zorder=4)
            dose_center_legend = \
                [Line2D([0], [0], color='r', ls='None', marker='+', \
                        label='Isodose Center')]
            ax.add_artist(ax.legend(handles=dose_center_legend, \
                          bbox_to_anchor=(1.25, 0.35), loc='upper left', \
                          borderaxespad=0.))
        
    return fig
        
def legend_dose_storage(folder_path):
    dose_storage_file = pydicom.dcmread(finder_rt_dose_file(folder_path))
    dose_plan_file = pydicom.dcmread(finder_rt_plan_file(folder_path))
    dose_color = [[0.55, 0, 1], [0, 0, 1], [0, 0.5, 1], [0, 1, 0],
                  [1, 1, 0], [1, 0.65, 0], [1, 0, 0], [0.55, 0, 0]]
    dose_level = [0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 1]
    dose_prescription = dose_plan_file[0x300a0010][0].TargetPrescriptionDose
    dose_value = [item*dose_prescription for item in dose_level]
    dose_name = []
    for item in range(0, len(dose_level)):
        dose_name.append(str(round(dose_level[item]*100, 2))+' % / ' +\
                      str(round(dose_value[item], 3))+' '+\
                      str(dose_storage_file.DoseUnits))  
    dose_storage_legend = {'dose color':dose_color, 'dose level':dose_level,\
                            'dose value':dose_value, 'dose name': dose_name}
    return dose_storage_legend

def legend_structure(structure_file):
    structure_total_number = len(structure_file.ROIContourSequence)
    structure_name = []
    edge_color = []
    face_color = []
    for item in range(0, structure_total_number):
        structure_name.append(structure_file.\
                              StructureSetROISequence[item].ROIName)
        edge_color.append(RGBNormalize(structure_file.\
                          ROIContourSequence[item].ROIDisplayColor))
        face_color.append(edge_color[item]+[0])
    structure_legend = {'structure name' : structure_name,\
                        'edge color' : edge_color,\
                        'face color' : face_color}
    return structure_legend
    
    
def RGBNormalize(rgb):
    """Normalize RGB values from [0, 255] to [0, 1]

        In the dcm files, the RGB values are [0, 255], 
        but in matplotlib are [0, 1]

        Args:
            rgb: The rgb list. [r, g, b]

        return:
            Normalized RGB values.
            [r_0to1, g_0to1, b_0to1]

    """
    normalized_rgb = []
    for item in range(0, 3):
        normalized_rgb.append(rgb[item]/255)
    return normalized_rgb     

def remove_2d_black_dot(pixel_array):
    for i in range(1, pixel_array.shape[0]-1):
        for j in range(1, pixel_array.shape[1]-1):
            if pixel_array[i][j] == 0:   
                neighbors = [pixel_array[i-1][j-1], pixel_array[i-1][j],\
                             pixel_array[i-1][j+1],\
                             pixel_array[i][j-11], pixel_array[i][j+1],\
                             pixel_array[i+1][j-1], pixel_array[i+1][j],\
                             pixel_array[i+1][j+1]]
                if neighbors.count(0) <= 3:
                    neighbors_nozero = []
                    for k in range(0, len(neighbors)):
                        if neighbors[k] != 0:
                            neighbors_nozero.append(neighbors[k])
                    pixel_array[i][j] = np.average(neighbors_nozero)
    return pixel_array
    
    

 
    

    
    
