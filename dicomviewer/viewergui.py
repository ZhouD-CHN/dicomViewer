#!/usr/bin/env pythonw
# -*- coding: utf-8 -*-
# Copyright (c) 2019 Zhou Dejun
# This file is the GUI for the dicom-rt file viewer prototype.
# The layout of this GUI was generated with the help of the software \
# wxFormBuilder.
# This file is available at https://github.com/McZhD/dicomViewer

import wx
import wx.xrc
from matplotlib.backends.backend_wxagg \
    import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg \
    import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.figure import Figure
import dicomrotate

class DicomImg(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent=parent, size=wx.Size(1000, 900))

        self.fig = Figure()
        self.canvas = FigureCanvas(self, -1, self.fig)
        self.NavigationToolbar = NavigationToolbar(self.canvas)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.NavigationToolbar, proportion=0, border=2, \
                       flag=wx.ALL | wx.EXPAND)
        self.sizer.Add(self.canvas, 1, wx.EXPAND)
        self.SetSizer(self.sizer)


    def showimage(self, folder_path, file_index, \
                  show_structure=False, show_isodose=False):
        self.fig.clf() 
        size_pixel = tuple(self.GetSize())
        figure =  dicomrotate.plot_axial_plane(folder_path, file_index,\
                showstructure = show_structure, showisodose = show_isodose)
        figure.set_figwidth(size_pixel[0]/100)
        figure.set_figheight(size_pixel[1]/100)
        self.canvas.figure = figure
        self.canvas.draw()       
        self.canvas.Refresh()

class ControlPanel(wx.Panel):
    def __init__(self, parent, top):
        wx.Panel.__init__(self, parent=parent, size=wx.Size(1000, 100))

        self.image = top

        topSizer = wx.FlexGridSizer(3, 1, 0, 0)
        topSizer.SetFlexibleDirection(wx.BOTH)
        topSizer.SetNonFlexibleGrowMode(wx.FLEX_GROWMODE_SPECIFIED)

        sizer1 = wx.FlexGridSizer(0, 6, 0, 0)
        sizer1.SetFlexibleDirection(wx.BOTH)
        sizer1.SetNonFlexibleGrowMode(wx.FLEX_GROWMODE_SPECIFIED)

        self.openfile = wx.Button(self, wx.ID_ANY, u"Open Patient",\
                                  wx.DefaultPosition, wx.DefaultSize, 0)
        sizer1.Add(self.openfile, 0, wx.ALL, 5)

        self.patientdir = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString,\
                                      wx.DefaultPosition, wx.Size(400, -1), 0)
        sizer1.Add(self.patientdir, 0, wx.ALL, 5)

        self.curslice = wx.StaticText(self, wx.ID_ANY, u"Current Slice", \
                                      wx.DefaultPosition, wx.DefaultSize, 0)
        self.curslice.Wrap(-1)

        sizer1.Add(self.curslice, 0, wx.ALL, 5)

        self.slicno = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString, \
                                  wx.DefaultPosition, wx.DefaultSize, 0)
        sizer1.Add(self.slicno, 0, wx.ALL, 5)

        self.totslice = wx.StaticText(self, wx.ID_ANY, u"Total Slice", \
                                      wx.DefaultPosition, wx.DefaultSize, 0)
        self.totslice.Wrap(-1)

        sizer1.Add(self.totslice, 0, wx.ALL, 5)

        self.totslicno = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString, \
                                     wx.DefaultPosition, wx.DefaultSize, 0)
        sizer1.Add(self.totslicno, 0, wx.ALL, 5)

        topSizer.Add(sizer1, 1, wx.EXPAND, 5)

        sizer2 = wx.FlexGridSizer(0, 5, 0, 0)
        sizer2.SetFlexibleDirection(wx.BOTH)
        sizer2.SetNonFlexibleGrowMode(wx.FLEX_GROWMODE_SPECIFIED)

        self.previous = wx.Button(self, wx.ID_ANY, u"Previous Slice", \
                                  wx.DefaultPosition, wx.DefaultSize, 0)
        sizer2.Add(self.previous, 0, wx.ALL, 5)

        self.next = wx.Button(self, wx.ID_ANY, u"Next Slice", \
                              wx.DefaultPosition, wx.DefaultSize, 0)
        sizer2.Add(self.next, 0, wx.ALL, 5)

        self.jumpto = wx.StaticText(self, wx.ID_ANY, u"Go to Slice:", \
                                    wx.DefaultPosition, wx.DefaultSize, 0)
        self.jumpto.Wrap(-1)

        sizer2.Add(self.jumpto, 0, wx.ALL, 5)

        self.jumpslicno = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString, \
                                      wx.DefaultPosition, wx.DefaultSize, 0)
        sizer2.Add(self.jumpslicno, 0, wx.ALL, 5)

        self.jump = wx.Button(self, wx.ID_ANY, u"Go", wx.DefaultPosition, \
                              wx.DefaultSize, 0)
        sizer2.Add(self.jump, 0, wx.ALL, 5)

        topSizer.Add(sizer2, 1, wx.EXPAND, 5)

        sizer3 = wx.FlexGridSizer(0, 3, 0, 0)
        sizer3.SetFlexibleDirection(wx.BOTH)
        sizer3.SetNonFlexibleGrowMode(wx.FLEX_GROWMODE_SPECIFIED)

        self.showstru = wx.CheckBox(self, wx.ID_ANY, u"Show Structure", \
                                    wx.DefaultPosition, wx.DefaultSize, 0)
        sizer3.Add(self.showstru, 0, wx.ALL, 5)

        self.showdose = wx.CheckBox(self, wx.ID_ANY, u"Show Isodose", \
                                    wx.DefaultPosition, wx.DefaultSize, 0)
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
        dlg = wx.DirDialog(self, u"Choose Patient Folder", \
                           style=wx.DD_DEFAULT_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            datapath = dlg.GetPath()
            self.totalno = 1
            self.imageno = 1
            self.slicno.ChangeValue(str(self.imageno))
            self.totalno = len(dicomrotate.finder_ct_image_file(datapath))
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
        self.image.showimage(self.patientdir.GetValue(), self.imageno,
                        self.showstru.IsChecked(), self.showdose.IsChecked())
        self.slicno.ChangeValue(self.jumpslicno.GetValue())
        self.imageno = int(self.jumpslicno.GetValue())

    def reload(self, event):
        self.image.showimage(self.patientdir.GetValue(), self.imageno,
                        self.showstru.IsChecked(), self.showdose.IsChecked())


class Main(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, parent=None, \
                    title="Dicom-RT File Viewer Prototype", size=(1000, 1000))

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

