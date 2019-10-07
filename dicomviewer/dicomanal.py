#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 Zhou Dejun
# This file contains functions to analyze '.dcm'.
# This file is available at https://github.com/McZhD/dicomViewer

import matplotlib.pyplot as plt
import numpy as np
import pydicom
import os


def RTStructureFileFinder(folder_path):
    """RT Structure Set Storage file Finder

        Locate the RT Structure Set Storage File in the patient folder.

        Args:
            folder_path: The path to the patient folder.

        Returns:
            If there is no structure file in the folder, returns a string. Else return
            the path to the RT Structure Set Storage file.
    """
    file_list = []
    structure_file_path = 'No Structure File'
    for fpath, dirs, fs in os.walk(folder_path):
        for f in fs:
            filename = os.path.join(fpath, f)
            if filename.endswith('.dcm'):
                file_list.append(filename)
    for item in file_list:
        if pydicom.dcmread(item).SOPClassUID.name == 'RT Structure Set Storage':
            structure_file_path = item

    return structure_file_path


def RTDoseFileFinder(folder_path):
    """RT Dose Storage file Finder

        Locate the RT Structure Set Storage File in the patient folder.

        Args:
            folder_path: The path to the patient folder.

        Returns:
            If there is no dose storage file in the folder, returns a string. Else return
            the path to the RT Dose Storage file.
    """
    file_list = []
    dose_file_path = 'No Dose File'
    for fpath, dirs, fs in os.walk(folder_path):
        for f in fs:
            filename = os.path.join(fpath, f)
            if filename.endswith('.dcm'):
                file_list.append(filename)
    for item in file_list:
        if pydicom.dcmread(item).SOPClassUID.name == 'RT Dose Storage':
            dose_file_path = item

    return dose_file_path


def RTPlanFileFinder(folder_path):
    """RT Plan Storage file Finder

        Locate the RT Structure Set Storage File in the patient folder.

        Args:
            folder_path: The path to the patient folder.

        Returns:
            If there is no dose plan file in the folder, returns a string. Else return
            the path to the RT Plan Storage file.
    """
    file_list = []
    plan_file_path = 'No Plan File'
    for fpath, dirs, fs in os.walk(folder_path):
        for f in fs:
            filename = os.path.join(fpath, f)
            if filename.endswith('.dcm'):
                file_list.append(filename)
    for item in file_list:
        if pydicom.dcmread(item).SOPClassUID.name == 'RT Plan Storage':
            plan_file_path = item

    return plan_file_path


def CTImageFileFinder(folder_path):
    """CT Image Storage files Finder

        Locate the CT Image Storage files in the patient folder.

        Args:
            folder_path: The path to the patient folder.

        Returns:
            Return a list of the path to the CT Image Storage files.
    """
    file_list = []
    ct_image_file_list = []
    for fpath, dirs, fs in os.walk(folder_path):
        for f in fs:
            filename = os.path.join(fpath, f)
            if filename.endswith('.dcm'):
                file_list.append(filename)
    file_list.sort()
    for item in file_list:
        if pydicom.dcmread(item).SOPClassUID.name == 'CT Image Storage':
            ct_image_file_list.append(item)

    return ct_image_file_list


def AxialSlicePixelArray(folder_path, slice_number):
    """Axial Slice Pixel Array

        Extract the Axial Slice Pixel Array of the designated slice.

        Args:
            folder_path: The path to the patient folder.
            slice_number: The No. of the CT slice in the return of def CTImageFileFinder, start from 1.

        Returns:
            Return a the XY plane Pixel Array of the designated CT slice and the patient position of the slice.
            ([number of rows],[number of columns]), [x, y, z]
    """
    if slice_number > len(CTImageFileFinder(folder_path)):
        return 'Slice Number Out Of Range'
    slice_data = pydicom.dcmread(CTImageFileFinder(folder_path)[slice_number-1])
    slice_pixel_array = slice_data.pixel_array
    slice_position = slice_data.ImagePositionPatient
    return slice_pixel_array, slice_position


def AxialSliceStructureInformation(folder_path, slice_number):
    """Axial Slice Structure Information

        Extract the Structure information the designated axial slice.

        Args:
            folder_path: The path to the patient folder.
            slice_number: The No. of the CT slice in the return of def CTImageFileFinder, start from 1.

        Returns:
            List of names of the structures, List of colors of the structures,
            List of contour coordinates of the structures
           [name1, name2, ...], [[r_name1, g_name1, b_name1,...],
           [[[[list of name1_contour1_x],[list of name1_contour1_y], [list of name1_contour1_z], ...], ...], ...]
    """
    structure_file_path = RTStructureFileFinder(folder_path)
    slice_position = pydicom.dcmread(CTImageFileFinder(folder_path)[slice_number - 1]).ImagePositionPatient
    pixel_spacing = pydicom.dcmread(CTImageFileFinder(folder_path)[slice_number - 1]).PixelSpacing

    if structure_file_path == 'No Structure File':
        print(structure_file_path)
        return structure_file_path

    if slice_number > len(CTImageFileFinder(folder_path)):
        return 'Slice Number Out Of Range'

    structure_data = pydicom.dcmread(structure_file_path)
    structure_total_number = len(structure_data.ROIContourSequence)
    structure_name = []
    structure_color = []
    structure_contour = [[] for i in range(structure_total_number)]

    for item in range(0, structure_total_number):
        structure_name.append(structure_data.StructureSetROISequence[item].ROIName)
        structure_color.append(RGBNormalize(structure_data.ROIContourSequence[item].ROIDisplayColor))
    for item in range(0, structure_total_number):
        contour_sequence = structure_data.ROIContourSequence[item].ContourSequence
        for ctr in range(0, len(contour_sequence)):
            if slice_position[2] == contour_sequence[ctr].ContourData[2]:
                temp = contour_sequence[ctr].ContourData
                tempx = []
                tempy = []
                tempz = []
                for i in range(0, int(len(temp)/3)):
                    tempx.append((temp[i * 3 + 0]-slice_position[0])/pixel_spacing[0])
                    tempy.append((temp[i * 3 + 1]-slice_position[1])/pixel_spacing[1])
                    tempz.append(temp[i * 3 + 2])
                structure_contour[item].append([tempx, tempy, tempz])

    return structure_name, structure_color, structure_contour


def AxialSliceDoseInformation(folder_path, slice_number):
    """Axial Slice Dose Information

        Extract the Dose Distribution information the designated axial slice.

        Args:
            folder_path: The path to the patient folder.
            slice_number: The No. of the CT slice in the return of def CTImageFileFinder, start from 1.

        Returns:
            Levels of dose, Colors of the levels, dose distribution, Labels of the dose.
            [level1, level2, ...], [color1, color2, ...], ([number of rows],[number of columns]), [label1, label2, ...]
    """
    dose_file_path = RTDoseFileFinder(folder_path)
    if dose_file_path == 'No Dose File':
        print(dose_file_path)
        return dose_file_path
    dose_data = pydicom.dcmread(dose_file_path)
    slice_data = pydicom.dcmread(CTImageFileFinder(folder_path)[slice_number-1])
    dose_grid_scale = dose_data.DoseGridScaling
    dose_max = np.max(dose_data.pixel_array)
    dose_unit = dose_data.DoseUnits
    dose_position = dose_data.ImagePositionPatient
    slice_thick = dose_data.SliceThickness
    dose_row = int(dose_data.Rows)
    dose_column = int(dose_data.Columns)
    dose_frame = int(dose_data.NumberOfFrames)
    dose_pixel_spacing = dose_data.PixelSpacing
    slice_position = slice_data.ImagePositionPatient
    slice_row = slice_data.Rows
    slice_column = slice_data.Columns

    dose_level = [0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 1]
    dose_value = [item*dose_max for item in dose_level]
    dose_name = []
    for item in range(0, len(dose_level)):
        dose_name.append(str(round(dose_level[item]*100, 2))+' % / ' +
                         str(round(dose_max * dose_grid_scale, 3))+' '+dose_unit)
    dose_color = [[0.55, 0, 1], [0, 0, 1], [0, 0.5, 1], [0, 1, 0],
                 [1, 1, 0], [1, 0.65, 0], [1, 0, 0], [0.55, 0, 0]]
    # dosectr = [[[], [], []] for i in range(len(doseperc))]
    dose_distribution = np.zeros([slice_row, slice_column])

    if slice_position[2] < dose_position[2] or slice_position[2] >= dose_position[2] + slice_thick*dose_frame:
        return dose_value, dose_color, dose_distribution, dose_name

    dose_pixel_array = dose_data.pixel_array[int((slice_position[2]-dose_position[2])/slice_thick)]
    row_start = int((dose_position[0]-slice_position[0])/dose_pixel_spacing[0])
    column_start = int((dose_position[1]-slice_position[1])/dose_pixel_spacing[1])
    for row in range(0, dose_row):
        for col in range(0, dose_column):
            dose_distribution[column_start+row][row_start+col] = dose_pixel_array[row][col]

    return dose_value, dose_color, dose_distribution, dose_name


def DrawAxialSlice(folder_path, slice_number, structure_show=False, dose_show=False):
    """Draw Axial Slice

        Draw and save the designated Axial slice.

        Args:
            folder_path: The path to the patient folder.
            slice_number: The No. of the CT slice in the return of def imgxyfile, start from 1.
            structure_show: Whether draw the contours of structure on the XY slice.
            dose_show: Whether draw the contours of dose distribution on the XY slice.

    """
    slice_pixel_array = AxialSlicePixelArray(folder_path, slice_number)[0]
    plt.figure(figsize=(11, 6))
    plt.title('Slice Position Z = ' + str(AxialSlicePixelArray(folder_path, slice_number)[1][2]))
    plt.imshow(slice_pixel_array, cmap=plt.bone(), zorder=1)

    if structure_show is True and AxialSliceStructureInformation(folder_path, slice_number) != 'No Structure File':
        structure_information = AxialSliceStructureInformation(folder_path, slice_number)
        for structure_name in range(0, len(structure_information[0])):
            color = structure_information[1][structure_name]
            alpha = (color[0], color[1], color[2], np.linspace(0.75, 0, len(structure_information[0]))[structure_name])
            for contour in range(0, len(structure_information[2][structure_name])):
                plt.fill(structure_information[2][structure_name][contour][0],
                         structure_information[2][structure_name][contour][1],
                         edgecolor=color, facecolor=alpha, label=structure_information[0][structure_name], zorder=2)
            handles, labels = plt.gca().get_legend_handles_labels()
            handle_str, label_str = [], []
            for handle, label in zip(handles, labels):
                if label not in label_str:
                    handle_str.append(handle)
                    label_str.append(label)
            structure_legend = plt.legend(handle_str, label_str, bbox_to_anchor=(-0.1, 1),
                                          loc='upper right', borderaxespad=0.)
            plt.gca().add_artist(structure_legend)
    if dose_show is True and AxialSliceDoseInformation(folder_path, slice_number) != 'No Dose File':
        dose_information = AxialSliceDoseInformation(folder_path, slice_number)
        isodose_lines = []
        cs = plt.contour(dose_information[2], dose_information[0], colors=dose_information[1], zorder=3)
        isodose_lines.extend(cs.collections)
        if len(isodose_lines) != 1:
            dose_legend = plt.legend(isodose_lines, dose_information[3],
                                     bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
            plt.gca().add_artist(dose_legend)
    plt.savefig(CTImageFileFinder(folder_path)[slice_number-1][0:-4]+'.png', dpi=100)
    plt.close()


def RGBNormalize(rgb):
    """Normalize RGB values from [0, 255] to [0, 1]

        In the dcm files, the RGB values are [0, 255], but in matplotlib are [0, 1]

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





