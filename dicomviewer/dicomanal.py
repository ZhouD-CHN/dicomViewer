#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 Zhou Dejun
# This file contains functions to analyze '.dcm'.
# This file is available at https://github.com/McZhD/dicomViewer

import matplotlib.pyplot as plt
import numpy as np
import pydicom
import os


def strufile(datadir):
    """RT Structure Set Storage file Finder

        Locate the RT Structure Set Storage File in the patient folder.

        Args:
            datadir: The path to the patient folder.

        Returns:
            If there is no structure file in the folder, returns a string. Else return
            the path to the RT Structure Set Storage file.
    """
    data = []
    strufilepath = 'No Structure File'
    for fpath, dirs, fs in os.walk(datadir):
        for f in fs:
            filename = os.path.join(fpath, f)
            if filename.endswith('.dcm'):
                data.append(filename)
    for item in data:
        if pydicom.dcmread(item).SOPClassUID.name == 'RT Structure Set Storage':
            strufilepath = item

    return strufilepath


def dosefile(datadir):
    """RT Dose Storage file Finder

        Locate the RT Structure Set Storage File in the patient folder.

        Args:
            datadir: The path to the patient folder.

        Returns:
            If there is no dose storage file in the folder, returns a string. Else return
            the path to the RT Dose Storage file.
    """
    data = []
    dosefilepath = 'No Dose File'
    for fpath, dirs, fs in os.walk(datadir):
        for f in fs:
            filename = os.path.join(fpath, f)
            if filename.endswith('.dcm'):
                data.append(filename)
    for item in data:
        if pydicom.dcmread(item).SOPClassUID.name == 'RT Dose Storage':
            dosefilepath = item

    return dosefilepath


def planfile(datadir):
    """RT Plan Storage file Finder

        Locate the RT Structure Set Storage File in the patient folder.

        Args:
            datadir: The path to the patient folder.

        Returns:
            If there is no dose plan file in the folder, returns a string. Else return
            the path to the RT Plan Storage file.
    """
    data = []
    planfilepath = 'No Plan File'
    for fpath, dirs, fs in os.walk(datadir):
        for f in fs:
            filename = os.path.join(fpath, f)
            if filename.endswith('.dcm'):
                data.append(filename)
    for item in data:
        if pydicom.dcmread(item).SOPClassUID.name == 'RT Plan Storage':
            planfilepath = item

    return planfilepath


def imgxyfile(datadir):
    """CT Image Storage files Finder

        Locate the CT Image Storage files in the patient folder.

        Args:
            datadir: The path to the patient folder.

        Returns:
            Return a list of the path to the CT Image Storage files.
    """
    data = []
    imgxylist = []
    for fpath, dirs, fs in os.walk(datadir):
        for f in fs:
            filename = os.path.join(fpath, f)
            if filename.endswith('.dcm'):
                data.append(filename)
    data.sort()
    for item in data:
        if pydicom.dcmread(item).SOPClassUID.name == 'CT Image Storage':
            imgxylist.append(item)

    return imgxylist


def imgxy(datadir, no):
    """Image XY Pixel Array

        Extract the XY plane Pixel Array of the designated CT slice.

        Args:
            datadir: The path to the patient folder.
            no: The No. of the CT slice in the return of def imgxyfile, start from 1.

        Returns:
            Return a the XY plane Pixel Array of the designated CT slice and the patient position of the slice.
            ([number of rows],[number of columns]), [x, y, z]
    """
    if no > len(imgxyfile(datadir)):
        return 'Slice Number Out Of Range'
    imgset = pydicom.dcmread(imgxyfile(datadir)[no-1])
    imgxyct = imgset.pixel_array
    imgpos = imgset.ImagePositionPatient
    return imgxyct, imgpos


def imgxystru(datadir, no):
    """Image XY Structure information

        Extract the Structure information the designated CT slice.

        Args:
            datadir: The path to the patient folder.
            no: The No. of the CT slice in the return of def imgxyfile, start from 1.

        Returns:
            List of names of the structures, List of colors of the structures,
            List of contour coordinates of the structures
           [name1, name2, ...], [[r_name1, g_name1, b_name1,...],
           [[[[list of name1_contour1_x],[list of name1_contour1_y], [list of name1_contour1_z], ...], ...], ...]
    """
    strupath = strufile(datadir)
    slicpos = pydicom.dcmread(imgxyfile(datadir)[no - 1]).ImagePositionPatient
    pixspc = pydicom.dcmread(imgxyfile(datadir)[no - 1]).PixelSpacing

    if strupath == 'No Structure File':
        print(strupath)
        return strupath

    if no > len(imgxyfile(datadir)):
        return 'Slice Number Out Of Range'

    struset = pydicom.dcmread(strupath)
    strutotno = len(struset.ROIContourSequence)
    struname = []
    strucolor = []
    structr = [[] for i in range(strutotno)]

    for item in range(0, strutotno):
        struname.append(struset.StructureSetROISequence[item].ROIName)
        strucolor.append(rgbnorm(struset.ROIContourSequence[item].ROIDisplayColor))
    for item in range(0, strutotno):
        ctrseq = struset.ROIContourSequence[item].ContourSequence
        for ctr in range(0, len(ctrseq)):
            if slicpos[2] == ctrseq[ctr].ContourData[2]:
                temp = ctrseq[ctr].ContourData
                tempx = []
                tempy = []
                tempz = []
                for i in range(0, int(len(temp)/3)):
                    tempx.append((temp[i * 3 + 0]-slicpos[0])/pixspc[0])
                    tempy.append((temp[i * 3 + 1]-slicpos[1])/pixspc[1])
                    tempz.append(temp[i * 3 + 2])
                structr[item].append([tempx, tempy, tempz])

    return struname, strucolor, structr


def imgxydose(datadir, no):
    """Image XY Dose Distribution information

        Extract the Dose Distribution information the designated CT slice.

        Args:
            datadir: The path to the patient folder.
            no: The No. of the CT slice in the return of def imgxyfile, start from 1.

        Returns:
            Levels of dose, Colors of the levels, dose distribution, Labels of the dose.
            [level1, level2, ...], [color1, color2, ...], ([number of rows],[number of columns]), [label1, label2, ...]
    """
    dosepath = dosefile(datadir)
    if dosepath == 'No Dose File':
        print(dosepath)
        return dosepath
    doseset = pydicom.dcmread(dosepath)
    slicset = pydicom.dcmread(imgxyfile(datadir)[no-1])
    dosegridscal = doseset.DoseGridScaling
    dosemax = np.max(doseset.pixel_array)
    doseunit = doseset.DoseUnits
    imgpos = doseset.ImagePositionPatient
    slicthic = doseset.SliceThickness
    slicpos = slicset.ImagePositionPatient
    slicrow = slicset.Rows
    sliccol = slicset.Columns
    imgrow = int(doseset.Rows)
    imgcol = int(doseset.Columns)
    imgfrm = int(doseset.NumberOfFrames)
    pixspc = doseset.PixelSpacing

    doseperc = [0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 1]
    doseval = [x*dosemax for x in doseperc]
    dosename = []
    for item in range(0, len(doseperc)):
        dosename.append(str(round(doseperc[item]*100, 2))+' % / '+str(round(dosemax * dosegridscal, 3))+' '+doseunit)
    dosecolor = [[0.55, 0, 1], [0, 0, 1], [0, 0.5, 1], [0, 1, 0],
                 [1, 1, 0], [1, 0.65, 0], [1, 0, 0], [0.55, 0, 0]]
    # dosectr = [[[], [], []] for i in range(len(doseperc))]
    doseimg = np.zeros([slicrow, sliccol])

    if slicpos[2] < imgpos[2] or slicpos[2] >= imgpos[2] + slicthic*imgfrm:
        return doseval, dosecolor, doseimg, dosename

    imgpix = doseset.pixel_array[int((slicpos[2]-imgpos[2])/slicthic)]
    rowstart = int((imgpos[0]-slicpos[0])/pixspc[0])
    colstart = int((imgpos[1]-slicpos[1])/pixspc[1])
    for row in range(0, imgrow):
        for col in range(0, imgcol):
            doseimg[colstart+row][rowstart+col] = imgpix[row][col]

    return doseval, dosecolor, doseimg, dosename


def drawimgxy(datadir, no, stru=False, dose=False):
    """Draw Image XY

        Draw and save the designated XY slice.

        Args:
            datadir: The path to the patient folder.
            no: The No. of the CT slice in the return of def imgxyfile, start from 1.
            stru: Whether draw the contours of structure on the XY slice.
            dose: Whether draw the contours of dose distribution on the XY slice.

    """
    slice = imgxy(datadir, no)[0]
    plt.figure(figsize=(11, 6))
    plt.title('Slice Position Z = ' + str(imgxy(datadir, no)[1][2]))
    plt.imshow(slice, cmap=plt.bone(), zorder=1)

    if stru is True and imgxystru(datadir, no) != 'No Structure File':
        structr = imgxystru(datadir, no)
        for struname in range(0, len(structr[0])):
            color = structr[1][struname]
            alpha = (color[0], color[1], color[2], np.linspace(0.75, 0, len(structr[0]))[struname])
            for ctr in range(0, len(structr[2][struname])):
                plt.fill(structr[2][struname][ctr][0], structr[2][struname][ctr][1],
                         edgecolor=color, facecolor=alpha, label=structr[0][struname], zorder=2)
            handles, labels = plt.gca().get_legend_handles_labels()
            handle_str, label_str = [], []
            for handle, label in zip(handles, labels):
                if label not in label_str:
                    handle_str.append(handle)
                    label_str.append(label)
            strulgd = plt.legend(handle_str, label_str, bbox_to_anchor=(-0.1, 1), loc='upper right', borderaxespad=0.)
            plt.gca().add_artist(strulgd)
    if dose is True and imgxydose(datadir, no) != 'No Dose File':
        dosectr = imgxydose(datadir, no)
        lines = []
        cs = plt.contour(dosectr[2], dosectr[0], colors=dosectr[1], zorder=3)
        lines.extend(cs.collections)
        if len(lines) != 1:
            doselgd = plt.legend(lines, dosectr[3], bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
            plt.gca().add_artist(doselgd)
    plt.savefig(imgxyfile(datadir)[no-1][0:-4]+'.png', dpi=100)
    plt.close()


def rgbnorm(rgb):
    """Normalize RGB values from [0, 255] to [0, 1]

        In the dcm files, the RGB values are [0, 255], but in matplotlib are [0, 1]

        Args:
            rgb: The rgb list. [r, g, b]

        return:
            Normalized RGB values.
            [r_0to1, g_0to1, b_0to1]

    """
    norm = []
    for item in range(0, 3):
        norm.append(rgb[item]/255)
    return norm





