#!/usr/bin/env pythonw
# -*- coding: utf-8 -*-
# Copyright (c) 2019 Zhou Dejun
# This file is the GUI for the dicom-rt file viewer prototype.
# The layout of this GUI was generated with the help of the software wxFormBuilder.
# This file is available at https://github.com/McZhD/dicomViewer

import wx
import wx.xrc
import numpy as np
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import dicomanal
#import dicomviewer.dicomanal as dicomanal


class DicomImg(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent=parent, size=wx.Size(1000, 700))

        self.fig = Figure()
        self.image = self.fig.add_subplot
        self.canvas = FigureCanvas(self, -1, self.fig)
        self.NavigationToolbar = NavigationToolbar(self.canvas)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.NavigationToolbar, proportion=0, border=2, flag=wx.ALL | wx.EXPAND)
        self.sizer.Add(self.canvas, 1, wx.EXPAND)
        self.SetSizer(self.sizer)

    def showimage(self, datadir, no, stru=False, dose=False):
        self.fig.clf()
        slice = dicomanal.imgxy(datadir, no)[0]
        title = 'Slice Position Z = ' + str(dicomanal.imgxy(datadir, no)[1][2])

        self.image().set_title(title)
        self.image().imshow(slice, cmap=plt.bone(), zorder=1)

        if stru is True and dicomanal.imgxystru(datadir, no) != 'No Structure File':
            structr = dicomanal.imgxystru(datadir, no)
            for struname in range(0, len(structr[0])):
                color = structr[1][struname]
                alpha = (color[0], color[1], color[2], np.linspace(0.75, 0, len(structr[0]))[struname])
                for ctr in range(0, len(structr[2][struname])):
                    self.image().fill(structr[2][struname][ctr][0], structr[2][struname][ctr][1],
                                      edgecolor=color, facecolor=alpha, label=structr[0][struname], zorder=2)
            handles, labels = self.image().get_legend_handles_labels()
            handle_str, label_str = [], []
            for handle, label in zip(handles, labels):
                if label not in label_str:
                    handle_str.append(handle)
                    label_str.append(label)
            strulgd = self.image().legend(handle_str, label_str,
                                          bbox_to_anchor=(-0.02, 1), loc='upper right', borderaxespad=0.)
            self.image().add_artist(strulgd)
        if dose is True and dicomanal.imgxydose(datadir, no) != 'No Dose File':
            dosectr = dicomanal.imgxydose(datadir, no)
            cs = self.image().contour(dosectr[2], dosectr[0], colors=dosectr[1], zorder=3)
            lines = []
            lines.extend(cs.collections)
            if len(lines) != 1:
                doselgd = self.image().legend(lines, dosectr[3],
                                              bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
                self.image().add_artist(doselgd)
        self.canvas.draw()


class ControlPanel(wx.Panel):
    def __init__(self, parent, top):
        wx.Panel.__init__(self, parent=parent, size=wx.Size(1000, 100))

        self.image = top

        # self.totalno = 1
        # self.imageno = 1

        topSizer = wx.FlexGridSizer(3, 1, 0, 0)
        topSizer.SetFlexibleDirection(wx.BOTH)
        topSizer.SetNonFlexibleGrowMode(wx.FLEX_GROWMODE_SPECIFIED)

        sizer1 = wx.FlexGridSizer(0, 6, 0, 0)
        sizer1.SetFlexibleDirection(wx.BOTH)
        sizer1.SetNonFlexibleGrowMode(wx.FLEX_GROWMODE_SPECIFIED)

        self.openfile = wx.Button(self, wx.ID_ANY, u"Open Patient", wx.DefaultPosition, wx.DefaultSize, 0)
        sizer1.Add(self.openfile, 0, wx.ALL, 5)

        self.patientdir = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size(400, -1), 0)
        sizer1.Add(self.patientdir, 0, wx.ALL, 5)

        self.curslice = wx.StaticText(self, wx.ID_ANY, u"Current Slice", wx.DefaultPosition, wx.DefaultSize, 0)
        self.curslice.Wrap(-1)

        sizer1.Add(self.curslice, 0, wx.ALL, 5)

        self.slicno = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0)
        sizer1.Add(self.slicno, 0, wx.ALL, 5)

        self.totslice = wx.StaticText(self, wx.ID_ANY, u"Total Slice", wx.DefaultPosition, wx.DefaultSize, 0)
        self.totslice.Wrap(-1)

        sizer1.Add(self.totslice, 0, wx.ALL, 5)

        self.totslicno = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0)
        sizer1.Add(self.totslicno, 0, wx.ALL, 5)

        topSizer.Add(sizer1, 1, wx.EXPAND, 5)

        sizer2 = wx.FlexGridSizer(0, 5, 0, 0)
        sizer2.SetFlexibleDirection(wx.BOTH)
        sizer2.SetNonFlexibleGrowMode(wx.FLEX_GROWMODE_SPECIFIED)

        self.previous = wx.Button(self, wx.ID_ANY, u"Previous Slice", wx.DefaultPosition, wx.DefaultSize, 0)
        sizer2.Add(self.previous, 0, wx.ALL, 5)

        self.next = wx.Button(self, wx.ID_ANY, u"Next Slice", wx.DefaultPosition, wx.DefaultSize, 0)
        sizer2.Add(self.next, 0, wx.ALL, 5)

        self.jumpto = wx.StaticText(self, wx.ID_ANY, u"Go to Slice:", wx.DefaultPosition, wx.DefaultSize, 0)
        self.jumpto.Wrap(-1)

        sizer2.Add(self.jumpto, 0, wx.ALL, 5)

        self.jumpslicno = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0)
        sizer2.Add(self.jumpslicno, 0, wx.ALL, 5)

        self.jump = wx.Button(self, wx.ID_ANY, u"Go", wx.DefaultPosition, wx.DefaultSize, 0)
        sizer2.Add(self.jump, 0, wx.ALL, 5)

        topSizer.Add(sizer2, 1, wx.EXPAND, 5)

        sizer3 = wx.FlexGridSizer(0, 2, 0, 0)
        sizer3.SetFlexibleDirection(wx.BOTH)
        sizer3.SetNonFlexibleGrowMode(wx.FLEX_GROWMODE_SPECIFIED)

        self.showstru = wx.CheckBox(self, wx.ID_ANY, u"Show Structure", wx.DefaultPosition, wx.DefaultSize, 0)
        sizer3.Add(self.showstru, 0, wx.ALL, 5)

        self.showdose = wx.CheckBox(self, wx.ID_ANY, u"Show Isodose", wx.DefaultPosition, wx.DefaultSize, 0)
        sizer3.Add(self.showdose, 0, wx.ALL, 5)

        topSizer.Add(sizer3, 1, wx.EXPAND, 5)

        self.openfile.Bind(wx.EVT_BUTTON, self.readpwd)
        self.next.Bind(wx.EVT_BUTTON, self.nextslic)
        self.previous.Bind(wx.EVT_BUTTON, self.preslic)
        self.jump.Bind(wx.EVT_BUTTON, self.goslic)
        self.showstru.Bind(wx.EVT_CHECKBOX, self.reload)
        self.showdose.Bind(wx.EVT_CHECKBOX, self.reload)

        self.SetSizer(topSizer)
        self.Layout()

    def readpwd(self, event):
        dlg = wx.DirDialog(self, u"Choose Patient Folder", style=wx.DD_DEFAULT_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            datapath = dlg.GetPath()
            self.totalno = 1
            self.imageno = 1
            self.slicno.ChangeValue(str(self.imageno))
            self.totalno = len(dicomanal.imgxyfile(datapath))
            self.totslicno.ChangeValue(str(self.totalno))
            self.patientdir.ChangeValue(datapath)
        dlg.Destroy()
        self.image.showimage(self.patientdir.GetValue(), self.imageno,
                             self.showstru.IsChecked(), self.showdose.IsChecked())

    def nextslic(self, event):
        self.imageno = int(self.imageno) + 1
        if self.imageno >= self.totalno:
            self.imageno = int(self.imageno) - self.totalno
        self.image.showimage(self.patientdir.GetValue(), self.imageno,
                             self.showstru.IsChecked(), self.showdose.IsChecked())
        self.slicno.ChangeValue(str(self.imageno))

    def preslic(self, event):
        self.imageno = int(self.imageno) - 1
        if self.imageno <= 0:
            self.imageno = int(self.imageno) + self.totalno
        self.image.showimage(self.patientdir.GetValue(), self.imageno,
                             self.showstru.IsChecked(), self.showdose.IsChecked())
        self.slicno.ChangeValue(str(self.imageno))

    def goslic(self, event):
        # self.imageno = int(self.jumpslicno.GetValue())
        self.image.showimage(self.patientdir.GetValue(), int(self.jumpslicno.GetValue()),
                             self.showstru.IsChecked(), self.showdose.IsChecked())
        self.slicno.ChangeValue(self.jumpslicno.GetValue())
        self.imageno = int(self.jumpslicno.GetValue())

    def reload(self, event):
        self.image.showimage(self.patientdir.GetValue(), self.imageno,
                             self.showstru.IsChecked(), self.showdose.IsChecked())


class Main(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, parent=None, title="Dicom-RT File Viewer Prototype", size=(1000, 800))

        splitter = wx.SplitterWindow(self)
        top = DicomImg(splitter)
        bottom = ControlPanel(splitter, top)
        splitter.SplitHorizontally(top, bottom)
        splitter.SetMinimumPaneSize(100)


if __name__ == "__main__":
    app = wx.App()
    frame = Main()
    frame.Show()
    app.MainLoop()

