#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import time
import re
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5 import QtCore, QtGui, QtWidgets
from scipy.odr import ODR, Model, RealData
from array import array
from pandas import DataFrame, concat

matplotlib.rcParams['agg.path.chunksize'] = 100000

class mplCustomizedToolbar(NavigationToolbar):
  """
  Modify the default matplotlib toolbar excluding some buttons that implicitly do the same that interface is designed to.
  """
  toolitems = [
      t for t in NavigationToolbar.toolitems
      if t[0] in ('Home', 'Pan', 'Zoom', 'Save')
  ]

  def __init__(self, *args, **kwargs):
    super(mplCustomizedToolbar, self).__init__(*args, **kwargs)
    self.layout().takeAt(4)


class MplCanvas(FigureCanvas):

  def __init__(self, parent=None, width=8.5, height=5.5, dpi=72):
    self.fig = Figure(figsize=(width, height), dpi=dpi)
    self.axes = self.fig.add_axes([0.075, 0.075, 0.85, 0.85],
                                  facecolor='w',
                                  label='plot')
    super(MplCanvas, self).__init__(self.fig)


class graphMainWindow(QtWidgets.QMainWindow):
  typeGraphMenu = ['line', 'scatter', 'histogram', 'autocorrelation']

  [
      grMenuItems, eijMenuItems, dataSet, labels, histogramTitles, plotTitles,
      dataTitles, xData, yData, nrData, viewAllCoord
  ] = [[] for _ in range(11)]

  extension = ''
  filename = ''
  integralMarker = False
  integralIndex = 0
  molDimAdjust = False
  horizontalGrLine = False

  Xlabel = ''
  Ylabel = ''
  GrLabel = ''
  EijLabel = ''

  changeX = False
  changeY = False
  changeHist = False

  # Information and values related to the last graph to be plotted.
  canvasInfo = {
      'type': '',
      'title': '',
      'data': [],
      'user parameters': '',
      'user data': [],
      'user boundaries': []
  }

  def __init__(self, parent=None):
    super().__init__(parent)
    #QtWidgets.QShortcut(QtWidgets.QKeySequence('Shift+Return'), self, self.plot)
    #self.setWindowIcon(QtWidgets.QIcon('icone.png'))
    #icone.png -> https://www.dropbox.com/s/t4kyvv3uhpb3a4a/icone.png
    self.setWindowTitle('DiceWin')
    self.upperMenuBar()
    self.status = QtWidgets.QStatusBar()
    self.setStatusBar(self.status)
    self.layoutElements()

  def layoutElements(self):
    """
    (None) -> None
    Set up the widgets that will compose the main window.
    """
    ## graph options widgets
    self.eijMenulLabel = QtWidgets.QLabel('Output list:')
    self.eijMenulLabel.hide()
    self.eijMenu = QtWidgets.QComboBox()
    self.eijMenu.addItems(self.eijMenuItems)
    self.eijMenu.setStyleSheet('QComboBox {combobox-popup: 0;}')
    self.eijMenu.setMaxVisibleItems(10)
    self.eijMenu.currentIndexChanged.connect(self.changeSimulationOutput)
    self.eijMenu.hide()

    self.dataXMenu = QtWidgets.QComboBox()
    self.dataXMenu.addItems(self.labels)
    self.dataXMenu.currentIndexChanged.connect(self.changeXData)

    self.dataYMenu = QtWidgets.QComboBox()
    self.dataYMenu.addItems(self.labels)
    self.dataYMenu.currentIndexChanged.connect(self.changeYData)

    self.intervalXMin = QtWidgets.QLineEdit()
    self.intervalXMax = QtWidgets.QLineEdit()
    self.intervalYMin = QtWidgets.QLineEdit()
    self.intervalYMax = QtWidgets.QLineEdit()

    self.plotXIntervalMin = QtWidgets.QLineEdit()
    self.plotXIntervalMax = QtWidgets.QLineEdit()
    self.plotYIntervalMin = QtWidgets.QLineEdit()
    self.plotYIntervalMax = QtWidgets.QLineEdit()

    ## histogram options widgets
    self.checkGaussian = QtWidgets.QCheckBox('Plot Gaussian')
    self.histIndexMin = QtWidgets.QDoubleSpinBox()
    self.histIndexMin.setDecimals(0)
    self.histIndexMin.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                    QtWidgets.QSizePolicy.Fixed)

    self.histIndexMax = QtWidgets.QDoubleSpinBox()
    self.histIndexMax.setDecimals(0)
    self.histIndexMax.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                    QtWidgets.QSizePolicy.Fixed)

    self.binValue = QtWidgets.QDoubleSpinBox()
    self.binValue.setDecimals(0)
    self.binValue.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                QtWidgets.QSizePolicy.Fixed)

    ## autocorrelation options widgets
    self.largestT = QtWidgets.QLineEdit()

    self.expMenu = QtWidgets.QComboBox()
    self.expMenu.addItems([
        'A1*exp(-t/B1)', 'A1*exp(-t/B1)+A2*exp(-t/B2)',
        'A1*exp(-t/B1)+A2*exp(-t/B2)+A3*exp(-t/B3)'
    ])
    self.expMenu.currentIndexChanged.connect(self.changeFitParameters)

    self.A1 = QtWidgets.QLineEdit()
    self.A1.setPlaceholderText('A1')
    self.A1.setEnabled(True)

    self.A2 = QtWidgets.QLineEdit()
    self.A2.setPlaceholderText('A2')
    self.A2.setEnabled(False)

    self.A3 = QtWidgets.QLineEdit()
    self.A3.setPlaceholderText('A3')
    self.A3.setEnabled(False)

    self.B1 = QtWidgets.QLineEdit()
    self.B1.setPlaceholderText('B1')
    self.B1.setEnabled(True)

    self.B2 = QtWidgets.QLineEdit()
    self.B2.setPlaceholderText('B2')
    self.B2.setEnabled(False)

    self.B3 = QtWidgets.QLineEdit()
    self.B3.setPlaceholderText('B3')
    self.B3.setEnabled(False)

    self.buttonApply = QtWidgets.QPushButton('Apply')
    self.buttonApply.clicked.connect(self.applyUserParameters)

    self.largestTLabel = QtWidgets.QLabel('')

    ## rdf widgets
    self.grMenu = QtWidgets.QComboBox()
    self.grMenu.addItems(self.grMenuItems)
    self.grMenu.setStyleSheet('QComboBox {combobox-popup: 0;}')
    self.grMenu.setMaxVisibleItems(10)
    self.grMenu.currentIndexChanged.connect(self.changeXData)
    self.grMenu.currentIndexChanged.connect(self.changeYData)

    self.molDimA = QtWidgets.QLineEdit()
    self.molDimB = QtWidgets.QLineEdit()
    self.molDimC = QtWidgets.QLineEdit()

    self.buttonApplyRDF = QtWidgets.QPushButton('Apply')
    self.buttonApplyRDF.clicked.connect(self.applyMolDim)

    ## ueff widgets
    self.temperature = QtWidgets.QLineEdit()
    self.temperature.setText('298.0')

    ## canvas widgets
    self.graphxyFrame = QtWidgets.QWidget()
    self.canvas = MplCanvas()
    self.cid = self.canvas.fig.canvas.mpl_connect('button_press_event',
                                                  self.clickIntegral)
    # self.mpl_toolbar = mplCustomizedToolbar(self.canvas, self.graphxyFrame)
    self.mpl_toolbar = NavigationToolbar(self.canvas, self.graphxyFrame)

    self.checkOverplot = QtWidgets.QCheckBox('Overplot')

    self.buttonRedefinePlotInterval = QtWidgets.QPushButton('Apply')
    self.buttonRedefinePlotInterval.clicked.connect(self.applyCanvasBoundaries)

    self.buttonViewAll = QtWidgets.QPushButton('View all')
    self.buttonViewAll.clicked.connect(self.viewAll)

    ## output widgets
    self.oa1 = QtWidgets.QLabel('A1=')
    self.oa2 = QtWidgets.QLabel('A2=')
    self.oa3 = QtWidgets.QLabel('A3=')
    self.ob1 = QtWidgets.QLabel('B1=')
    self.ob2 = QtWidgets.QLabel('B2=')
    self.ob3 = QtWidgets.QLabel('B3=')

    self.fittingLabel = QtWidgets.QLabel('fitting<br/>parameters')

    self.integralLabel = QtWidgets.QLabel('')

    self.meanLabel = QtWidgets.QLabel('')

    self.stdDeviationLabel = QtWidgets.QLabel('')

    self.typeGraph = QtWidgets.QComboBox()
    self.typeGraph.addItems(self.typeGraphMenu)
    self.typeGraph.currentIndexChanged.connect(self.changeDefaultBoundaries)

    self.buttonPlot = QtWidgets.QPushButton('Show')
    self.buttonPlot.clicked.connect(self.plot)

    self.buttonSave = QtWidgets.QPushButton('Save')
    self.buttonSave.clicked.connect(self.savePlotData)

    self.intervalXMin.returnPressed.connect(self.chcondX)
    self.intervalXMax.returnPressed.connect(self.chcondX)
    self.intervalYMin.returnPressed.connect(self.chcondY)
    self.intervalYMax.returnPressed.connect(self.chcondY)

    self.plotXIntervalMin.returnPressed.connect(self.applyCanvasBoundaries)
    self.plotXIntervalMax.returnPressed.connect(self.applyCanvasBoundaries)
    self.plotYIntervalMin.returnPressed.connect(self.applyCanvasBoundaries)
    self.plotYIntervalMax.returnPressed.connect(self.applyCanvasBoundaries)

    self.A1.returnPressed.connect(self.applyUserParameters)
    self.A2.returnPressed.connect(self.applyUserParameters)
    self.A3.returnPressed.connect(self.applyUserParameters)
    self.B1.returnPressed.connect(self.applyUserParameters)
    self.B2.returnPressed.connect(self.applyUserParameters)
    self.B3.returnPressed.connect(self.applyUserParameters)
    self.largestT.returnPressed.connect(self.applyUserParameters)

    self.molDimA.returnPressed.connect(self.applyMolDim)
    self.molDimA.setPlaceholderText('A')
    self.molDimB.returnPressed.connect(self.applyMolDim)
    self.molDimB.setPlaceholderText('B')
    self.molDimC.returnPressed.connect(self.applyMolDim)
    self.molDimC.setPlaceholderText('C')

    #
    AllEntries = [
        self.intervalXMin, self.intervalXMax, self.intervalYMin,
        self.intervalYMax, self.plotXIntervalMin, self.plotXIntervalMax,
        self.plotYIntervalMin, self.plotYIntervalMax, self.largestT, self.A1,
        self.A2, self.A3, self.B1, self.B2, self.B3, self.molDimA, self.molDimB,
        self.molDimC, self.temperature
    ]
    for e in AllEntries:
      e.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                      QtWidgets.QSizePolicy.Minimum)

    self.horizontalLine = QtWidgets.QFrame()
    self.horizontalLine.setFrameShape(QtWidgets.QFrame.HLine)
    self.horizontalLine.setFrameShadow(QtWidgets.QFrame.Sunken)
    self.blankWidget = QtWidgets.QWidget()

    ## graph options frame
    self.gridGraphEntries = QtWidgets.QGridLayout()
    self.gridGraphEntries.addWidget(self.eijMenulLabel, 0, 0, 1, 1,
                                    QtCore.Qt.AlignHCenter)
    self.gridGraphEntries.addWidget(self.eijMenu, 0, 1, 1, 4)
    self.gridGraphEntries.addWidget(QtWidgets.QLabel('<small>x-axis</small>'),
                                    1, 1, 1, 1, QtCore.Qt.AlignHCenter)
    self.gridGraphEntries.addWidget(QtWidgets.QLabel('<small>y-axis</small>'),
                                    1, 4, 1, 1, QtCore.Qt.AlignHCenter)
    self.gridGraphEntries.addWidget(QtWidgets.QLabel('Data'), 2, 0, 1, 1,
                                    QtCore.Qt.AlignHCenter)
    self.gridGraphEntries.addWidget(self.dataXMenu, 2, 1, 1, 1)
    self.gridGraphEntries.addWidget(QtWidgets.QLabel('Min. value:'), 3, 0, 1, 1,
                                    QtCore.Qt.AlignRight)
    self.gridGraphEntries.addWidget(QtWidgets.QLabel('Max. value:'), 4, 0, 1, 1,
                                    QtCore.Qt.AlignRight)
    self.gridGraphEntries.addWidget(self.intervalXMin, 3, 1, 1, 1)
    self.gridGraphEntries.addWidget(self.intervalXMax, 4, 1, 1, 1)
    self.gridGraphEntries.addWidget(QtWidgets.QLabel('<b>Data<b/>'), 2, 3, 1, 1,
                                    QtCore.Qt.AlignHCenter)
    self.gridGraphEntries.addWidget(self.dataYMenu, 2, 4, 1, 1)
    self.gridGraphEntries.addWidget(QtWidgets.QLabel('Min. value:'), 3, 3, 1, 1,
                                    QtCore.Qt.AlignRight)
    self.gridGraphEntries.addWidget(QtWidgets.QLabel('Max. value:'), 4, 3, 1, 1,
                                    QtCore.Qt.AlignRight)
    self.gridGraphEntries.addWidget(self.intervalYMin, 3, 4, 1, 1)
    self.gridGraphEntries.addWidget(self.intervalYMax, 4, 4, 1, 1)

    self.widgetButtonsEntries = QtWidgets.QGroupBox('Line and Scatter options')
    self.widgetButtonsEntries.setLayout(self.gridGraphEntries)

    ## histogram options frame
    self.gridHistogramEntries = QtWidgets.QGridLayout()
    self.gridHistogramEntries.addWidget(self.checkGaussian, 0, 0, 1, 2,
                                        QtCore.Qt.AlignCenter)
    self.gridHistogramEntries.addWidget(QtWidgets.QLabel('Num. of bins:'), 1, 0,
                                        1, 1, QtCore.Qt.AlignRight)
    self.gridHistogramEntries.addWidget(self.binValue, 1, 1, 1, 1,
                                        QtCore.Qt.AlignVCenter)
    self.gridHistogramEntries.addWidget(QtWidgets.QLabel('Min. index:'), 0, 2,
                                        1, 1, QtCore.Qt.AlignRight)
    self.gridHistogramEntries.addWidget(self.histIndexMin, 0, 3, 1, 1)
    self.gridHistogramEntries.addWidget(QtWidgets.QLabel('Max. index:'), 1, 2,
                                        1, 1, QtCore.Qt.AlignRight)
    self.gridHistogramEntries.addWidget(self.histIndexMax, 1, 3, 1, 1)

    self.widgetHistogramEntries = QtWidgets.QGroupBox('Histogram options')
    self.widgetHistogramEntries.setLayout(self.gridHistogramEntries)

    ## autocorrelation options frame
    self.gridAutoCorrEntries = QtWidgets.QGridLayout()
    self.gridAutoCorrEntries.addWidget(QtWidgets.QLabel('Fit type:'), 1, 0, 1,
                                       1, QtCore.Qt.AlignRight)
    self.gridAutoCorrEntries.addWidget(self.expMenu, 1, 1, 1, 4)
    self.gridAutoCorrEntries.addWidget(QtWidgets.QLabel('A:'), 2, 0, 1, 1,
                                       QtCore.Qt.AlignRight)
    self.gridAutoCorrEntries.addWidget(self.A1, 2, 1, 1, 1)
    self.gridAutoCorrEntries.addWidget(self.A2, 2, 2, 1, 1)
    self.gridAutoCorrEntries.addWidget(self.A3, 2, 3, 1, 1)
    self.gridAutoCorrEntries.addWidget(QtWidgets.QLabel('B:'), 3, 0, 1, 1,
                                       QtCore.Qt.AlignRight)
    self.gridAutoCorrEntries.addWidget(self.B1, 3, 1, 1, 1)
    self.gridAutoCorrEntries.addWidget(self.B2, 3, 2, 1, 1)
    self.gridAutoCorrEntries.addWidget(self.B3, 3, 3, 1, 1)
    self.gridAutoCorrEntries.addWidget(QtWidgets.QLabel('Largest t:'), 4, 0, 1,
                                       1, QtCore.Qt.AlignRight)
    self.gridAutoCorrEntries.addWidget(self.largestT, 4, 1, 1, 1)
    self.gridAutoCorrEntries.addWidget(self.largestTLabel, 4, 2, 1, 2)
    self.gridAutoCorrEntries.addWidget(self.buttonApply, 4, 3, 1, 1)

    self.widgetAutoCorrEntries = QtWidgets.QGroupBox('Autocorrelation options')
    self.widgetAutoCorrEntries.setLayout(self.gridAutoCorrEntries)

    ## rdf options frame
    self.gridGr = QtWidgets.QGridLayout()
    self.gridGr.addWidget(self.grMenu, 0, 0, 1, 3)
    self.gridGr.addWidget(
        QtWidgets.QLabel('<small>MOLECULAR DIMENSIONS</small>'), 1, 0, 1, 3,
        QtCore.Qt.AlignLeft)
    #self.gridGr.addWidget(QtWidgets.QLabel('A:'), 2, 0, 1, 1, QtCore.Qt.AlignRight)
    self.gridGr.addWidget(self.molDimA, 2, 0, 1, 1)
    #self.gridGr.addWidget(QtWidgets.QLabel('B:'), 2, 2, 1, 1, QtCore.Qt.AlignRight)
    self.gridGr.addWidget(self.molDimB, 2, 1, 1, 1)
    #self.gridGr.addWidget(QtWidgets.QLabel('C:'), 2, 4, 1, 1, QtCore.Qt.AlignRight)
    self.gridGr.addWidget(self.molDimC, 2, 2, 1, 1)
    self.gridGr.addWidget(self.buttonApplyRDF, 3, 2, 1, 1)
    self.gridGr.addWidget(QtWidgets.QLabel('<b>Ueff options</b>'), 4, 0, 1, 3)
    self.gridGr.addWidget(QtWidgets.QLabel('Temperature:'), 5, 0, 1, 1,
                          QtCore.Qt.AlignRight)
    self.gridGr.addWidget(self.temperature, 5, 1, 1, 1)
    self.gridGr.addWidget(QtWidgets.QLabel('K'), 5, 2, 1, 1)
    #self.gridGr.setColumnStretch(0, 2)
    #self.gridGr.setColumnStretch(1, 1)
    #self.gridGr.setColumnStretch(2, 2)

    self.widgetGr = QtWidgets.QGroupBox('RDF options')
    self.widgetGr.setLayout(self.gridGr)
    self.widgetGr.hide()

    ## Ueff options frame
    #self.gridUeff = QtWidgets.QUeffidLayout()
    #self.gridUeff.addWidget(QtWidgets.QLabel('Temperature'), 0, 0, 1, 1)
    #self.gridUeff.addWidget(self.temperature, 0, 1, 1, 1)
    #self.gridUeff.addWidget(QtWidgets.QLabel('K'), 0, 3, 1, 1)

    #self.widgetUeff = QtWidgets.QUeffoupBox('Ueff options')
    #self.widgetUeff.setLayout(self.gridUeff)
    #self.widgetUeff.hide()

    ## type of graph
    self.gridTypeOfGraph = QtWidgets.QGridLayout()
    self.gridTypeOfGraph.addWidget(self.typeGraph, 0, 0, 1, 1)
    self.gridTypeOfGraph.addWidget(self.buttonPlot, 0, 1, 1, 1)
    self.gridTypeOfGraph.addWidget(self.buttonSave, 0, 2, 1, 1)
    self.gridTypeOfGraph.setColumnStretch(0, 2)
    self.gridTypeOfGraph.setColumnStretch(1, 1)
    self.gridTypeOfGraph.setColumnStretch(2, 1)

    self.widgetTypeOfGraph = QtWidgets.QGroupBox('Graph options')
    self.widgetTypeOfGraph.setLayout(self.gridTypeOfGraph)

    ## output frame
    self.gridOutput = QtWidgets.QGridLayout()
    self.gridOutput.addWidget(self.meanLabel, 1, 0, 1, 2)
    self.gridOutput.addWidget(self.stdDeviationLabel, 1, 2, 1, 2)
    self.gridOutput.addWidget(self.integralLabel, 0, 0, 1, 1)
    self.gridOutput.addWidget(self.fittingLabel, 2, 0, 2, 1,
                              QtCore.Qt.AlignHCenter)
    self.gridOutput.addWidget(self.oa1, 2, 1, 1, 1)
    self.gridOutput.addWidget(self.oa2, 2, 2, 1, 1)
    self.gridOutput.addWidget(self.oa3, 2, 3, 1, 1)
    self.gridOutput.addWidget(self.ob1, 3, 1, 1, 1)
    self.gridOutput.addWidget(self.ob2, 3, 2, 1, 1)
    self.gridOutput.addWidget(self.ob3, 3, 3, 1, 1)

    self.widgetOutput = QtWidgets.QFrame()
    self.widgetOutput.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.widgetOutput.setLayout(self.gridOutput)

    ## graph intervals frame
    self.gridAxesIntervals = QtWidgets.QGridLayout()
    self.gridAxesIntervals.addWidget(QtWidgets.QLabel('Min. X:'), 0, 0, 1, 1,
                                     QtCore.Qt.AlignRight)
    self.gridAxesIntervals.addWidget(self.plotXIntervalMin, 0, 1, 1, 1)
    self.gridAxesIntervals.addWidget(QtWidgets.QLabel('Max. X:'), 0, 2, 1, 1,
                                     QtCore.Qt.AlignRight)
    self.gridAxesIntervals.addWidget(self.plotXIntervalMax, 0, 3, 1, 1)
    self.gridAxesIntervals.addWidget(QtWidgets.QLabel('Min. Y:'), 1, 0, 1, 1,
                                     QtCore.Qt.AlignRight)
    self.gridAxesIntervals.addWidget(self.plotYIntervalMin, 1, 1, 1, 1)
    self.gridAxesIntervals.addWidget(QtWidgets.QLabel('Max. Y:'), 1, 2, 1, 1,
                                     QtCore.Qt.AlignRight)
    self.gridAxesIntervals.addWidget(self.plotYIntervalMax, 1, 3, 1, 1)
    self.gridAxesIntervals.addWidget(self.buttonRedefinePlotInterval, 0, 4, 2,
                                     1)

    self.widgetAxesIntervals = QtWidgets.QGroupBox('Axes range options')
    self.widgetAxesIntervals.setLayout(self.gridAxesIntervals)

    ## canvas frame
    self.boxGraph = QtWidgets.QGridLayout()
    self.boxGraph.addWidget(self.mpl_toolbar, 0, 0, 1, 4)
    # self.boxGraph.addWidget(self.coordsText, 0, 2, 1, 1, QtCore.Qt.AlignRight)
    self.boxGraph.addWidget(self.horizontalLine, 1, 0, 1, 4)
    self.boxGraph.addWidget(self.canvas, 2, 0, 1, 4)
    self.boxGraph.addWidget(self.checkOverplot, 3, 0, 1, 1)
    self.boxGraph.addWidget(self.buttonViewAll, 3, 2, 1, 1,
                            QtCore.Qt.AlignRight)
    self.boxGraph.addWidget(self.widgetAxesIntervals, 6, 0, 1, 4)
    self.boxGraph.setColumnStretch(0, 2)
    self.boxGraph.setColumnStretch(1, 1)
    self.boxGraph.setColumnStretch(2, 2)

    self.frameGraph = QtWidgets.QFrame()
    self.frameGraph.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frameGraph.setLayout(self.boxGraph)
    self.frameGraph.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                  QtWidgets.QSizePolicy.Minimum)

    ##
    self.boxOuterFrame = QtWidgets.QGridLayout()
    self.boxOuterFrame.addWidget(self.widgetButtonsEntries, 0, 0, 1, 1)
    self.boxOuterFrame.addWidget(self.widgetHistogramEntries, 1, 0, 1, 1)
    self.boxOuterFrame.addWidget(self.widgetAutoCorrEntries, 2, 0, 1, 1)
    self.boxOuterFrame.addWidget(self.widgetGr, 2, 0, 1, 1)
    self.boxOuterFrame.addWidget(self.widgetTypeOfGraph, 3, 0, 1, 1)
    self.boxOuterFrame.addWidget(self.widgetOutput, 4, 0, 1, 1)
    self.boxOuterFrame.addWidget(self.frameGraph, 0, 1, 5, 1)

    self.mainFrame = QtWidgets.QWidget()
    self.mainFrame.setLayout(self.boxOuterFrame)
    self.mainFrame.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                 QtWidgets.QSizePolicy.Minimum)

    self.setCentralWidget(self.mainFrame)
    self.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                       QtWidgets.QSizePolicy.Minimum)

  def upperMenuBar(self):
    """
    (None) -> None
    Set up the widgets that will compose the upper menu.
    """
    upperMenuBar = self.menuBar()

    menuFile = upperMenuBar.addMenu('&File')

    openAction = QtWidgets.QAction('Open...', self)
    openAction.setShortcut('Ctrl+O')
    openAction.triggered.connect(self.selectDataFile)
    menuFile.addAction(openAction)
    menuFile.addSeparator()

    saveDataAs = QtWidgets.QAction('Save data', self)
    saveDataAs.setShortcut('Ctrl+S')
    saveDataAs.triggered.connect(self.savePlotData)
    menuFile.addAction(saveDataAs)
    menuFile.addSeparator()

    closeAction = QtWidgets.QAction('Close', self)
    closeAction.setShortcut('Ctrl+C')
    menuFile.addAction(closeAction)

    menuAbout = upperMenuBar.addMenu('&About')
    infoAction = QtWidgets.QAction('Info', self)
    infoAction.triggered.connect(self.aboutWindow)
    menuAbout.addAction(infoAction)

  def chcondX(self):
    self.changeX = True
    self.plot()

  def chcondY(self):
    self.changeY = True
    self.plot()

  def aboutWindow(self):
    """
    (None) -> None
    Create a pop-up window containing information about the interface.
    """
    info = QtWidgets.QMessageBox()
    info.setWindowTitle('About')
    info.setText('Application written in Python<br/><br/> \
                  Authors: Thiago Duarte, Emanuel Mancio and Kaline Coutinho<br/> \
                  Institution: Physics Institute, University of Sao Paulo<br/> \
                  Funding: CNPq-PIBIC and FAPESP<br/> \
                  Year: 2015 and 2020<br/><br/> \
                  Last Modifyd date: April 1, 2021')
    info.setIcon(1)
    info.exec_()

  def read_block(self, f, labels):
    data = {lab: array('f') for lab in labels}

    enum = list(enumerate(labels))

    new_block = False

    for line in f:
      line_values = line.split()
      try:
        for j, lab in enum:
          data[lab].append(float(line_values[j]))
      except:
        if line_values[0] == 'NMOVE':
          new_block = True
          break
        else:
          pass    # ignore line with data ****

    data = DataFrame(data)

    return data, new_block

  def readFile(self):
    """
    (None) -> List, List, List, List, String, String
    Create a window that allows the user to select the file which will be read. Return the data read and file name information.
    """
    openFile = QtWidgets.QFileDialog()
    openFile.setDirectory(os.getcwd())
    openFile.setFileMode(QtWidgets.QFileDialog.ExistingFile)
    openFile.setViewMode(0)
    selectedFileName = openFile.getOpenFileName()[0]

    if selectedFileName:
      labels, grItems, eijItems = [], [], []
      filename = str((selectedFileName.split('.'))[0])
      extension = str((selectedFileName.split('.'))[-1])

      # *.out files
      if (extension == 'out'):
        with open(selectedFileName, 'rt') as f:
          for line in f:
            if "MC steps" in line:
              sim_len = int(line.split()[-1])
            if 'NMOVE' in line:
              labels = line.rstrip().split()
              f.readline()
              break

          if sim_len >= 1000000:
            self.status.showMessage("Reading this file may take some minutes")
          data = {
              lab: np.zeros(sim_len, dtype=np.float32) for lab in labels[1:]
          }

          enum = list(enumerate(labels[1:]))

          i = -1
          sim_len -= 1

          for line in f:
            try:
              if line[-5] == "#":
                i += 1
                line_values = line.split()[1:-1]
                for j, lab in enum:
                  data[lab][i] = line_values[j]
              else:
                pass
            except:
              self.status.clearMessage()
              break

          if i != sim_len:
            self.status.showMessage("This simulation didn't finished")

        data = DataFrame(data)
        data["NMOVE"] = np.arange(1, sim_len + 2, dtype=np.uint32)
        data = data.reindex(columns=labels, copy=False)

      # *.dst files
      elif (extension == 'dst'):
        with open(selectedFileName, 'rt') as f:
          labels = f.readline().split()
          data = {lab: array('f') for lab in labels}
          for line in f:
            try:
              line_values = line.split()
              for i in range(len(line)):
                data[labels[i]].append(float(line_values[i]))
            except:
              pass    # if one value in the line is "****", the line won't be used

          labels = [lab for lab in labels if lab not in ('i', 'j')]
          data.pop('i')
          data.pop('j')
          data = DataFrame(data)

      # *.hbd files
      elif (extension == 'hbd'):
        with open(selectedFileName, 'rt') as f:
          criteria = f.readline()[2:]
          labels = (f.readline().split())[4:]

          data = {lab: array('f') for lab in labels}
          for line in f:
            try:
              line_values = [line[29:35],line[35:42],line[42:51],line[51:60],line[60:69],line[69:78],line[78:87],line[87:96],line[96:]]
              for i, key in enumerate(data):
                data[key].append(float(line_values[i]))
            except:
              pass    # if one value in the line is "****", the line won't be used

          data = DataFrame(data)

      # *.gr files
      elif extension == 'gr':
        with open(selectedFileName, 'rt') as f:
          f.readline()
          labels = ['r', 'G(r)', 'N(r)']
          data = {lab: array('f') for lab in labels}
          enum = list(enumerate(data))
          blocks = []
          for line in f:
            linelist = line.split()
            try:
              if linelist[1] == "RDF":
                grlabel = '{}{}({})-{}{}({})'.format(linelist[4], linelist[5],
                                                     linelist[9], linelist[11],
                                                     linelist[12], linelist[16])
                grItems.append(grlabel)
              else:
                for i, lab in enum:
                  data[lab].append(float(linelist[i]))
            except:
              blocks.append(DataFrame(data))
              data = {lab: array('f') for lab in labels}

          data = concat(blocks, axis=1, keys=grItems)

      # *.eij files
      elif (((extension)[0]) == 'e') and ((((extension)[-2:]) == 'ij') or ((
          (extension)[-2:]).isdigit() == True)):
        with open(selectedFileName, 'rt') as f:
          labels = f.readline().split()
          count = 1
          eijItems.append('Simulation output {}'.format(count))
          data, new_block = self.read_block(f, labels)

          blocks = [data]
          while new_block:
            data, new_block = self.read_block(f, labels)
            count += 1
            eijItems.append('Simulation output {}'.format(count))
            blocks.append(data)

          data = concat(blocks, axis=1, keys=[c for c in range(1, count + 1)])

      elif extension == "xvg":
        with open(selectedFileName, 'rt') as f:
          labels = ['x']
          for line in f:
            if '#' in line:
              pass
            elif '@' in line:
              if (re.search(r's\d', line) != None) and ("legend" in line):
                labels.append(re.search(r'"(.*?)"', line).group()[1:-1])
            else:
              initval = line.split()
              break

          if len(labels) == 1:
            val_qt = len(initval)
            for i in range(val_qt-1):
              labels.append("y{}".format(i))

          data = {lab: array('f') for lab in labels}
          enum = list(enumerate(labels))
          
          for i, lab in enum:
            data[lab].append(float(initval[i]))

          for line in f:
            values = line.split()
            if len(values) == len(labels): # pass empty or incomplete lines
              for i, lab in enum:
                data[lab].append(float(values[i]))

          data = DataFrame(data)

      # *.avr and generic files
      else:
        with open(selectedFileName, 'rt') as f:
          labels = gen_new_labels(f.readline().split())
          data = {lab: array('f') for lab in labels}
          for line in f:
            try:
              line_values = line.split()
              for i, key in enumerate(data):
                data[key].append(float(line_values[i]))
            except:
              pass    # if one value in the line is "****", the line won't be used

          data = DataFrame(data)

    try:
      return labels, data, grItems, eijItems, extension, filename
    except (UnboundLocalError):
      None

  def selectDataFile(self):
    """
    (None) -> None
    Fill the main lists with the data read from files and change layout according to the selected file extension.
    """
    try:
      self.labels, self.dataSet, self.grMenuItems, self.eijMenuItems, self.extension, self.filename = self.readFile(
      )

      if (self.extension == 'gr'):
        self.typeGraphMenu = ['line', 'scatter', 'histogram', 'ueff']
        self.default('gr')
        self.widgetGr.show()
        self.widgetAutoCorrEntries.hide()
        self.eijMenu.hide()
        self.eijMenulLabel.hide()
        self.fittingLabel.hide()

      elif (len(self.eijMenuItems) != 0):
        self.typeGraphMenu = ['line', 'scatter', 'histogram', 'autocorrelation']

        self.eijMenu.setCurrentIndex(0)
        self.default('normal')
        self.widgetAutoCorrEntries.show()
        self.widgetGr.hide()
        self.eijMenu.show()
        self.eijMenulLabel.show()
        self.fittingLabel.show()

      else:
        self.typeGraphMenu = ['line', 'scatter', 'histogram', 'autocorrelation']
        self.default('normal')
        self.widgetAutoCorrEntries.show()
        self.widgetGr.hide()
        self.eijMenu.hide()
        self.eijMenulLabel.hide()
        self.fittingLabel.show()

      self.oa1.show()
      self.oa2.show()
      self.oa3.show()
      self.ob1.show()
      self.ob2.show()
      self.ob3.show()

      # gmw = graphMainWindow()
      self.resize(self.minimumSizeHint())
    except TypeError as e:
      self.status.showMessage('ERROR: failed to open file.', 3456)

  def default(self, plotType):
    """
    (String) -> None
    Set up the default values and show the first graph according to the selected file extension.
    """

    self.layoutElements()

    if (plotType == 'normal'):
      self.dataXMenu.setEnabled(True)
      self.EijLabel = 1
      self.integralLabel.hide()

    elif (plotType == 'gr'):
      self.dataXMenu.setEnabled(False)
      self.integralLabel.show()

      self.grMenu.setCurrentIndex(1)

    self.dataXMenu.setCurrentIndex(0)
    self.dataYMenu.setCurrentIndex(1)

    xmin, xmax, ymin, ymax, ylen = self.updateData()
    self.intervalXMin.setText(str(xmin))
    self.intervalXMax.setText(str(xmax))
    self.intervalYMin.setText(str(ymin))
    self.intervalYMax.setText(str(ymax))

    self.histIndexMin.setMinimum(0)
    self.histIndexMin.setMaximum(ylen - 1)
    self.histIndexMin.setValue(0)
    self.histIndexMax.setMinimum(0)
    self.histIndexMax.setMaximum(ylen - 1)
    self.histIndexMax.setValue(ylen - 1)
    self.binValue.setMinimum(1)
    self.binValue.setMaximum(999)
    self.binValue.setValue(50)

    self.plotXIntervalMin.setText(str(xmin))
    self.plotXIntervalMax.setText(str(xmax))
    self.plotYIntervalMin.setText(str(ymin))
    self.plotYIntervalMax.setText(str(ymax))

    self.plot()

  def updateData(self):
    """
    (None) -> None
    Fill the menus with the data read from selected file.
    """
    try:
      if (self.extension == 'gr'):
        i = int(self.grMenu.currentIndex())
        k = int(self.dataYMenu.currentIndex())

        self.GrLabel = self.grMenuItems[i]
        self.Xlabel = 'r'
        self.Ylabel = self.labels[k]
        self.nrData = self.dataSet.loc[:, (self.GrLabel, 'N(r)')]
        xmin, xmax = self.dataSet.loc[:, (
            self.GrLabel,
            self.Xlabel)].min(), self.dataSet.loc[:, (self.GrLabel,
                                                      self.Xlabel)].max()
        ymin, ymax = self.dataSet.loc[:, (
            self.GrLabel,
            self.Ylabel)].min(), self.dataSet.loc[:, (self.GrLabel,
                                                      self.Ylabel)].max()
        ylen = len(self.dataSet.loc[:, (self.GrLabel, self.Ylabel)])

      elif (self.extension[0] == 'e'):
        i = int(self.eijMenu.currentIndex())
        j = int(self.dataXMenu.currentIndex())
        k = int(self.dataYMenu.currentIndex())

        self.EijLabel = i + 1
        self.Xlabel = self.labels[j]
        self.Ylabel = self.labels[k]
        xmin, xmax = self.dataSet.loc[:, (
            self.EijLabel,
            self.Xlabel)].min(), self.dataSet.loc[:, (self.EijLabel,
                                                      self.Xlabel)].max()
        ymin, ymax = self.dataSet.loc[:, (
            self.EijLabel,
            self.Ylabel)].min(), self.dataSet.loc[:, (self.EijLabel,
                                                      self.Ylabel)].max()
        ylen = len(self.dataSet.loc[:, (self.EijLabel, self.Ylabel)])

      else:
        j = int(self.dataXMenu.currentIndex())
        k = int(self.dataYMenu.currentIndex())
        self.Xlabel = self.labels[j]
        self.Ylabel = self.labels[k]
        xmin, xmax = self.dataSet.loc[:, self.Xlabel].min(
        ), self.dataSet.loc[:, self.Xlabel].max()
        ymin, ymax = self.dataSet.loc[:, self.Ylabel].min(
        ), self.dataSet.loc[:, self.Ylabel].max()
        ylen = len(self.dataSet.loc[:, self.Ylabel])

      return xmin, xmax, ymin, ymax, ylen
    except IndexError:
      self.status.showMessage(
          'ERROR: unable to load or update data. Please try to check if any file was already loaded.',
          3456)

  def dataProcessing(self):
    """
    (None) -> List
    Filter data from datasets according to values from 'Line and Scatter options'.
    """
    try:
      xMIN = float(self.intervalXMin.text())
      xMAX = float(self.intervalXMax.text())
      yMIN = float(self.intervalYMin.text())
      yMAX = float(self.intervalYMax.text())
      IDhMIN = int(self.histIndexMin.value())
      IDhMAX = int(self.histIndexMax.value())

      if (xMIN <= xMAX) and (yMIN <= yMAX):
        QtWidgets.QApplication.setOverrideCursor(
            QtGui.QCursor(QtCore.Qt.WaitCursor))

        if (self.extension == 'gr'):
          dataplot = self.dataSet.loc[:, self.GrLabel].copy()
          if (self.canvasInfo['type'] == 'histogram'):
            dataplot = dataplot.iloc[IDhMIN:IDhMAX, :]
          elif self.changeX or self.changeY:
            xData, yData = dataplot.loc[:,
                                        self.Xlabel], dataplot.loc[:,
                                                                   self.Ylabel]
            dataplot = dataplot.loc[(xData >= xMIN) & (xData <= xMAX) &
                                    (yData >= yMIN) & (yData <= yMAX)]
          elif self.canvasInfo['type'] == "ueff":
            dataplot[self.Ylabel] = dataplot[self.Ylabel].apply(self.ueff)

        elif self.extension[0] == 'e':
          dataplot = self.dataSet.loc[:, self.EijLabel]
          if self.canvasInfo["type"] == "autocorrelation":
            return dataplot
          elif (self.canvasInfo['type'] == 'histogram'):
            dataplot = dataplot.iloc[IDhMIN:IDhMAX, :]
          elif self.changeX or self.changeY:
            xData, yData = dataplot.loc[:,
                                        self.Xlabel], dataplot.loc[:,
                                                                   self.Ylabel]
            dataplot = dataplot.loc[(xData >= xMIN) & (xData <= xMAX) &
                                    (yData >= yMIN) & (yData <= yMAX)]

        else:
          dataplot = self.dataSet
          if self.canvasInfo["type"] == "autocorrelation":
            return dataplot
          elif (self.canvasInfo['type'] == 'histogram'):
            dataplot = dataplot.iloc[IDhMIN:IDhMAX, :]
          elif self.changeX or self.changeY:
            dataplot = self.dataSet.loc[:, [self.Ylabel, self.Xlabel]]
            xData, yData = dataplot.loc[:,
                                        self.Xlabel], self.dataSet.loc[:, self.
                                                                       Ylabel]
            dataplot = self.dataSet.loc[(xData >= xMIN) & (xData <= xMAX) &
                                        (yData >= yMIN) & (yData <= yMAX)]

        QtWidgets.QApplication.restoreOverrideCursor()
        return dataplot
      else:
        self.status.showMessage(
            'ERROR: invalid data interval range, Min. value > Max. value.',
            3456)
        return self.dataSet
    except (IndexError, ValueError):
      self.status.showMessage('ERROR: invalid data interval values.', 3456)
      return 1

  def plot(self):
    """
    (None) -> None
    Show graph according the user options.
    """
    try:
      graphHold = self.checkOverplot.isChecked()
      self.canvasInfo['type'] = str(self.typeGraph.currentText())
      self.integralMarker = False

      xid = self.dataXMenu.currentIndex()
      yid = self.dataYMenu.currentIndex()
      if (self.extension == 'gr'):
        gid = self.grMenu.currentIndex()
      if (self.extension[0] == 'e'):
        sid = self.eijMenu.currentIndex()

      dataplot = self.dataProcessing()
      data = []

      if (graphHold == False):
        self.plotTitles, self.histogramTitles = [], []
        self.canvas.axes.cla()

      if (self.extension
          == 'gr') and (self.canvasInfo['type'] != 'autocorrelation'):
        self.canvas.axes.set_position([0.075, 0.075, 0.7, 0.85])
      else:
        self.canvas.axes.set_position([0.075, 0.075, 0.85, 0.85])

      # WIDGETS RESPONSES
      if self.canvasInfo['type'] != 'histogram':
        self.meanLabel.hide()
        self.stdDeviationLabel.hide()
      if self.canvasInfo['type'] == 'autocorrelation':
        self.oa1.setText('A1=')
        self.oa2.setText('A2=')
        self.oa3.setText('A3=')
        self.ob1.setText('B1=')
        self.ob2.setText('B2=')
        self.ob3.setText('B3=')
        self.checkOverplot.hide()
      else:
        self.checkOverplot.show()

      # Line and Scatter
      if (self.canvasInfo['type']
          == 'line') or (self.canvasInfo['type']
                         == 'scatter') or (self.canvasInfo['type'] == 'ueff'):
        if (self.canvasInfo['type'] == 'line') or (self.canvasInfo['type']
                                                   == 'ueff'):
          self.canvas.axes.plot(self.Xlabel,
                                self.Ylabel,
                                data=dataplot,
                                linestyle='-')
        elif (self.canvasInfo['type'] == 'scatter'):
          self.canvas.axes.plot(self.Xlabel,
                                self.Ylabel,
                                data=dataplot,
                                linestyle='',
                                marker='.')
        self.horizontalUeffLine = False
        if (self.extension == 'gr'):
          self.canvasInfo['data'] = dataplot
          if (self.horizontalGrLine == False
              or graphHold == False and self.canvasInfo['type'] != 'ueff'):
            self.canvas.axes.axhline(y=1,
                                     c='k',
                                     linestyle='-',
                                     label='_nolegend_')
            self.horizontalGrLine = True
          if (self.horizontalUeffLine == False or graphHold == False):
            self.canvas.axes.axhline(y=0,
                                     c='k',
                                     linestyle='-',
                                     label='_nolegend_')
            self.horizontalUeffLine = True
          self.molDimAdjust = False
          self.plotTitles.append(r'${}$'.format(self.GrLabel))
          self.canvasInfo['title'] = '{}'.format(self.GrLabel)
          self.canvas.axes.legend(self.plotTitles,
                                  bbox_to_anchor=(1.05, 1),
                                  loc=2,
                                  ncol=1,
                                  mode='expand',
                                  borderaxespad=0.,
                                  frameon=False,
                                  numpoints=1,
                                  prop={'size': 10},
                                  handlelength=0.6,
                                  borderpad=-0.8)
        else:
          self.canvasInfo['data'] = dataplot
          xLabel = (self.Xlabel).replace('_', '-')
          yLabel = (self.Ylabel).replace('_', '-')
          self.plotTitles.append(r'${}\, \times\, {}$'.format(xLabel, yLabel))
          self.canvasInfo['title'] = '{}-{}'.format(xLabel, yLabel)
          self.canvas.axes.legend(self.plotTitles,
                                  bbox_to_anchor=(0., 1.052, 1., .9),
                                  loc=3,
                                  ncol=4,
                                  mode='expand',
                                  borderaxespad=0.,
                                  frameon=False,
                                  numpoints=1,
                                  prop={'size': 12},
                                  handlelength=0.6,
                                  borderpad=-0.8)

        xmin, xmax = self.canvas.axes.xaxis.get_data_interval()
        ymin, ymax = self.canvas.axes.yaxis.get_data_interval()
        self.setCanvasBoundaries(xmin, xmax, ymin, ymax)

      # Histogram
      elif (self.canvasInfo['type'] == 'histogram'):
        title = (self.Ylabel).replace('_', '-')
        self.histogramTitles.append(r'${}$'.format(title))

        try:
          numbins = int(self.binValue.value())
          if (numbins <= 0):
            numbins = 50
            self.status.showMessage(
                'ERROR: invalid number of bins. It was replaced by the default value.',
                3456)
        except (ValueError):
          numbins = 50

        L = len(dataplot.loc[:, self.Ylabel])
        _, bins = np.histogram(dataplot.loc[:, self.Ylabel], numbins)
        weight = [(1 / L) for _ in dataplot.loc[:, self.Ylabel]]

        n, bins, _ = self.canvas.axes.hist(dataplot.loc[:, self.Ylabel],
                                           numbins,
                                           weights=weight)

        if (self.extension == 'gr'):
          self.horizontalGrLine = False
          self.horizontalUeffLine = False
          self.canvas.axes.legend(self.histogramTitles,
                                  bbox_to_anchor=(1.05, 1),
                                  loc=2,
                                  ncol=1,
                                  mode='expand',
                                  borderaxespad=0.,
                                  frameon=False,
                                  numpoints=1,
                                  prop={'size': 10},
                                  handlelength=0.6,
                                  borderpad=-0.8)
        else:
          self.canvas.axes.legend(self.histogramTitles,
                                  bbox_to_anchor=(0.05, 1.055, 0.95, .9),
                                  loc=3,
                                  ncol=4,
                                  mode='expand',
                                  borderaxespad=0.,
                                  frameon=False,
                                  numpoints=1,
                                  prop={'size': 12},
                                  handlelength=0.6,
                                  borderpad=-0.8)

        binwidth = (bins[1] - bins[0])
        x = (min(bins) - binwidth)
        X = max(bins) + binwidth
        y = 0
        Y = max(n)

        decPrec = len((str(dataplot.loc[:,
                                        self.Ylabel].iloc[-1])).split('.')[-1])
        mean = float(dataplot.loc[:, self.Ylabel].mean())
        stdDev = float(dataplot.loc[:, self.Ylabel].std(ddof=0))

        points = np.linspace(float(bins[0]), float(bins[-1]), 250)

        if self.checkGaussian.isChecked():
          self.canvas.axes.plot(
              points,
              binwidth * ((1 / (stdDev * np.sqrt(2 * 3.14))) * np.exp(-0.5 * (
                  (points - mean) / stdDev)**2)),
              c='k',
              label='_nolegend_')

        self.meanLabel.show()
        self.stdDeviationLabel.show()
        self.meanLabel.setText('Mean: {}'.format(round(mean, decPrec)))
        self.stdDeviationLabel.setText('Standard deviation: {}'.format(
            round(stdDev, decPrec)))

        binCenter = [(bins[i - 1] + abs(bins[i] - bins[i - 1]) * 0.5)
                     for i in range(1, len(bins))]

        xmin, xmax = self.canvas.axes.xaxis.get_data_interval()
        ymin, ymax = self.canvas.axes.yaxis.get_data_interval()
        self.setCanvasBoundaries(xmin, xmax, ymin, ymax)
        self.canvasInfo['title'] = '{}'.format(title)
        self.canvasInfo['data'] = [n, binCenter]
        self.canvasInfo['user parameters'] = [mean, stdDev]
        self.canvasInfo['user data'] = [
            binwidth *
            ((1 /
              (stdDev * np.sqrt(2 * 3.14))) * np.exp(-0.5 *
                                                     ((i - mean) / stdDev)**2))
            for i in binCenter
        ]

      # Autocorrelation
      elif (self.canvasInfo['type'] == 'autocorrelation'):
        yLabel = (self.labels[yid]).replace('_', '-')
        self.plotTitles.append(r'$C(t) \times t\,\,\, ({})$'.format(yLabel))

        QtWidgets.QApplication.setOverrideCursor(
            QtGui.QCursor(QtCore.Qt.WaitCursor))
        x = np.arange(len(dataplot.loc[:, self.Ylabel]))
        y = self.estimated_autocorrelation(dataplot.loc[:,
                                                        self.Ylabel].to_numpy())
        QtWidgets.QApplication.restoreOverrideCursor()

        self.canvasInfo['title'] = '{}'.format(yLabel)
        self.canvasInfo['data'] = [x, y]

        self.canvas.axes.plot(self.canvasInfo['data'][0],
                              self.canvasInfo['data'][1],
                              linestyle='',
                              marker='.')
        self.setCanvasBoundaries(min(x), max(x), min(y), max(y))

        self.canvas.axes.legend(self.plotTitles,
                                bbox_to_anchor=(0.05, 1.055, 0.95, .9),
                                loc=3,
                                ncol=4,
                                mode='expand',
                                borderaxespad=0.,
                                frameon=False,
                                numpoints=1,
                                prop={'size': 12},
                                handlelength=0.6,
                                borderpad=-0.8)

      self.applyCanvasBoundaries()
      self.updatePlotLimits()

      self.canvas.draw()

    except (AttributeError) as e:
      # print(e)
      self.status.showMessage(
          'ERROR: failed to plot data. Please check the axes range values.',
          3456)
    except (IndexError) as e:
      # print(e)
      self.status.showMessage(
          'ERROR: failed to plot data. Please check data and axes range values.',
          3456)

  def applyUserParameters(self):
    """
    (None) -> None
    Get values from 'Autocorrelation options' and send to the fitting function.
    """
    if (self.canvasInfo['type'] == 'autocorrelation'):
      n = int(self.expMenu.currentIndex())
      largT = int(self.largestT.text())
      a1 = self.A1.text()
      a2 = self.A2.text()
      a3 = self.A3.text()
      b1 = self.B1.text()
      b2 = self.B2.text()
      b3 = self.B3.text()
      par = [a1, b1, a2, b2, a3, b3]
      self.bestFit(largT, par, n)
    else:
      self.status.showMessage(
          'ERROR: to change the guess it is necessary show an autocorrelation graph first.',
          3456)

  def bestFit(self, LT, P, i):
    """
    (Int, List, Int) -> None
    Fit two exponential curves the points generated in the autocorrelation function, one using coeficients inserted
    by the user and another using coeficients generated by a fitting function.
    """
    try:
      QtWidgets.QApplication.setOverrideCursor(
          QtGui.QCursor(QtCore.Qt.WaitCursor))
      if (i == 0):
        fit = Model(self.exponential1)
        guess = [float(k) for k in P[:2]]
      if (i == 1):
        fit = Model(self.exponential2)
        guess = [float(k) for k in P[:4]]
      if (i == 2):
        fit = Model(self.exponential3)
        guess = [float(k) for k in P]

      dataReal = RealData(self.canvasInfo['data'][0][:LT],
                          self.canvasInfo['data'][1][:LT])
      modelOdr = ODR(dataReal, fit, guess)
      modelOdr.set_job(fit_type=2)
      modelOutput = modelOdr.run()

      xnew = np.linspace(0, LT, len(self.canvasInfo['data'][0]))
      if (i == 0):
        yfit = self.exponential1(modelOutput.beta, xnew)
        yuser = self.exponential1(guess, xnew)
        initial = self.exponential1(guess, self.canvasInfo['data'][0][:LT])
        bestfit = self.exponential1(modelOutput.beta,
                                    self.canvasInfo['data'][0][:LT])
        self.oa1.setText('A1={}'.format(round(modelOutput.beta[0], 4)))
        self.ob1.setText('B1={}'.format(round(modelOutput.beta[1], 4)))
      if (i == 1):
        yfit = self.exponential2(modelOutput.beta, xnew)
        yuser = self.exponential2(guess, xnew)
        initial = self.exponential2(guess, self.canvasInfo['data'][0][:LT])
        bestfit = self.exponential2(modelOutput.beta,
                                    self.canvasInfo['data'][0][:LT])
        self.oa1.setText('A1={}'.format(round(modelOutput.beta[0], 4)))
        self.ob1.setText('B1={}'.format(round(modelOutput.beta[1], 4)))
        self.oa2.setText('A2={}'.format(round(modelOutput.beta[2], 4)))
        self.ob2.setText('B2={}'.format(round(modelOutput.beta[3], 4)))
      if (i == 2):
        yfit = self.exponential3(modelOutput.beta, xnew)
        yuser = self.exponential3(guess, xnew)
        initial = self.exponential3(guess, self.canvasInfo['data'][0][:LT])
        bestfit = self.exponential3(modelOutput.beta,
                                    self.canvasInfo['data'][0][:LT])
        self.oa1.setText('A1={}'.format(round(modelOutput.beta[0], 4)))
        self.ob1.setText('B1={}'.format(round(modelOutput.beta[1], 4)))
        self.oa2.setText('A2={}'.format(round(modelOutput.beta[2], 4)))
        self.ob2.setText('B2={}'.format(round(modelOutput.beta[3], 4)))
        self.oa3.setText('A3={}'.format(round(modelOutput.beta[4], 4)))
        self.ob3.setText('B3={}'.format(round(modelOutput.beta[5], 4)))

      # adjust boundaries
      self.canvas.axes.plot(self.canvasInfo['data'][0][:LT + 1],
                            self.canvasInfo['data'][1][:LT + 1],
                            linestyle='-')

      self.canvas.axes.clear()
      self.canvas.axes.plot(self.canvasInfo['data'][0],
                            self.canvasInfo['data'][1],
                            linestyle='',
                            marker='.',
                            label=r'$C(t)$')
      x, X = self.canvas.axes.get_xlim()
      y, Y = self.canvas.axes.get_ylim()
      self.canvas.axes.plot(xnew,
                            yuser,
                            linestyle='--',
                            c='k',
                            linewidth='1.5',
                            label=r'$guess$')
      self.canvas.axes.plot(xnew,
                            yfit,
                            linestyle='-',
                            c='r',
                            linewidth='1.5',
                            label=r'$fit$')

      self.canvas.axes.legend(bbox_to_anchor=(0., 1., 1., 1.),
                              loc=3,
                              ncol=3,
                              mode='expand',
                              borderaxespad=0.,
                              frameon=False,
                              numpoints=1,
                              prop={'size': 12},
                              handlelength=0.5)

      xmin, xmax = self.canvas.axes.xaxis.get_data_interval()
      ymin, ymax = self.canvas.axes.yaxis.get_data_interval()

      self.setCanvasBoundaries(xmin, LT, ymin, ymax)
      self.applyCanvasBoundaries()
      self.canvas.draw()

      self.canvasInfo['user data'] = [
          self.canvasInfo['data'][0][:LT], self.canvasInfo['data'][1][:LT],
          initial, bestfit
      ]
      self.canvasInfo['user parameters'] = [[
          i if i != '' else 'None' for i in P
      ], modelOutput.beta,
                                            self.expMenu.currentText(), LT]

      self.viewAllCoord = [x, X, y, Y]

      QtWidgets.QApplication.restoreOverrideCursor()

    except (IndexError) as e:
      self.status.showMessage(
          'ERROR: invalid fitting parameters. Please check your guess then try again.',
          3456)
      QtWidgets.QApplication.restoreOverrideCursor()
    except (ValueError) as e:
      self.status.showMessage(
          'ERROR: invalid fitting parameters. Please check your guess then try again.',
          3456)
      QtWidgets.QApplication.restoreOverrideCursor()

  def ueff(self, val):
    if val < 1e-20:
      return 0
    else:
      return (-1.985e-3 * float(self.temperature.text()) * np.log(val))

  def exponential1(self, C, t):
    return C[0] * np.exp(-t / C[1])

  def exponential2(self, C, t):
    return C[0] * np.exp(-t / C[1]) + C[2] * np.exp(-t / C[3])

  def exponential3(self, C, t):
    return C[0] * np.exp(-t / C[1]) + C[2] * np.exp(-t / C[3]) + C[4] * np.exp(
        -t / C[5])

  def estimated_autocorrelation(self, x):
    """
    (Array) -> Array
    Calculate the autocorrelation. Algorithm from: http://stackoverflow.com/q/14297012/190597
    """
    n = len(x)
    variance = x.var()
    x = x - x.mean()
    r = np.correlate(x, x, mode='full')[-n:]
    result = r / (variance * (np.arange(n, 0, -1)))
    return result

  def changeXData(self):
    """
    (None) -> None
    Update interval entries and datasets when 'data X' menu has the selected index changed.
    """
    try:
      if (self.extension == 'gr'):
        i = int(self.grMenu.currentIndex())
        self.GrLabel = self.grMenuItems[i]
        self.Xlabel = 'r'
        xmin, xmax = self.dataSet.loc[:, (
            self.GrLabel,
            self.Xlabel)].min(), self.dataSet.loc[:, (self.GrLabel,
                                                      self.Xlabel)].max()

      elif (self.extension[0] == 'e'):
        j = int(self.dataXMenu.currentIndex())
        self.Xlabel = self.labels[j]
        xmin, xmax = self.dataSet.loc[:, (
            self.EijLabel,
            self.Xlabel)].min(), self.dataSet.loc[:, (self.EijLabel,
                                                      self.Xlabel)].max()

      else:
        j = int(self.dataXMenu.currentIndex())
        self.Xlabel = self.labels[j]
        xmin, xmax = self.dataSet.loc[:, self.Xlabel].min(
        ), self.dataSet.loc[:, self.Xlabel].max()

      self.intervalXMin.setText(str(xmin))
      self.intervalXMax.setText(str(xmax))
      self.plotXIntervalMin.setText(str(xmin))
      self.plotXIntervalMax.setText(str(xmax))
    except (IndexError):
      self.status.showMessage(
          'ERROR: failed to automatic change x-axis values.', 3456)

  def changeYData(self):
    """
    (None) -> None
    Update interval entries and datasets when 'data Y' menu has the selected index changed. 
    """
    try:
      if (self.extension == 'gr'):
        i = int(self.grMenu.currentIndex())
        k = int(self.dataYMenu.currentIndex())
        self.GrLabel = self.grMenuItems[i]
        self.Ylabel = self.labels[k]
        self.nrData = self.dataSet.loc[:, (self.GrLabel, 'N(r)')]
        ymin, ymax = self.dataSet.loc[:, (
            self.GrLabel,
            self.Ylabel)].min(), self.dataSet.loc[:, (self.GrLabel,
                                                      self.Ylabel)].max()
        ylen = len(self.dataSet.loc[:, (self.GrLabel, self.Ylabel)])

      elif (self.extension[0] == 'e'):
        j = int(self.dataYMenu.currentIndex())
        self.Ylabel = self.labels[j]
        ymin, ymax = self.dataSet.loc[:, (
            self.EijLabel,
            self.Ylabel)].min(), self.dataSet.loc[:, (self.EijLabel,
                                                      self.Ylabel)].max()
        ylen = len(self.dataSet.loc[:, (self.EijLabel, self.Ylabel)])

      else:
        k = int(self.dataYMenu.currentIndex())
        self.Ylabel = self.labels[k]
        ymin, ymax = self.dataSet.loc[:, self.Ylabel].min(
        ), self.dataSet.loc[:, self.Ylabel].max()
        ylen = len(self.dataSet.loc[:, self.Ylabel])

      self.intervalYMin.setText(str(ymin))
      self.intervalYMax.setText(str(ymax))
      self.plotYIntervalMin.setText(str(ymin))
      self.plotYIntervalMax.setText(str(ymax))

      self.largestTLabel.setText('(from {})'.format(ylen))
      self.largestT.setText('{}'.format(int(ylen * 0.01)))
      self.histIndexMax.setMinimum(0)
      self.histIndexMax.setMaximum(ylen - 1)
      self.histIndexMin.setValue(0)
      self.histIndexMax.setValue(ylen - 1)

    except (IndexError):
      self.status.showMessage(
          'ERROR: failed to automatic change y-axis values.', 3456)

  def changeFitParameters(self):
    """
    (None) -> None
    Activate or deactivate fitting parameters entries according to selected function in the 'Fit type' menu.
    """
    index = self.expMenu.currentIndex()
    if index == 0:
      self.A1.setEnabled(True)
      self.A2.setEnabled(False)
      self.A3.setEnabled(False)
      self.B1.setEnabled(True)
      self.B2.setEnabled(False)
      self.B3.setEnabled(False)
    elif index == 1:
      self.A1.setEnabled(True)
      self.A2.setEnabled(True)
      self.A3.setEnabled(False)
      self.B1.setEnabled(True)
      self.B2.setEnabled(True)
      self.B3.setEnabled(False)
    elif index == 2:
      self.A1.setEnabled(True)
      self.A2.setEnabled(True)
      self.A3.setEnabled(True)
      self.B1.setEnabled(True)
      self.B2.setEnabled(True)
      self.B3.setEnabled(True)

  def changeDefaultBoundaries(self):
    """
    (None) -> None
    Set new values for 'Axes range options' when the type of graph is changed.
    """
    if (self.extension != ''):
      typeGraph = self.typeGraph.currentText()
      if (typeGraph == 'line') or (typeGraph == 'scatter'):
        xmin, xmax, ymin, ymax, _ = self.updateData()
        self.plotXIntervalMin.setText(str(xmin))
        self.plotXIntervalMax.setText(str(xmax))
        self.plotYIntervalMin.setText(str(ymin))
        self.plotYIntervalMax.setText(str(ymax))

  def updatePlotLimits(self):
    """
    (None) -> None
    Update the list 'viewAllCoord' which stores boundary values of all graphs on the canvas.
    """
    if self.checkOverplot.isChecked() == True:
      try:
        xmin, xmax = self.canvas.axes.get_xlim()
        ymin, ymax = self.canvas.axes.get_ylim()
        if xmin < self.viewAllCoord[0]:
          self.viewAllCoord[0] = xmin
        if xmax > self.viewAllCoord[1]:
          self.viewAllCoord[1] = xmax
        if ymin < self.viewAllCoord[2]:
          self.viewAllCoord[2] = ymin
        if ymax > self.viewAllCoord[3]:
          self.viewAllCoord[3] = ymax
      except (IndexError) as e:
        xmin, xmax = self.canvas.axes.get_xlim()
        ymin, ymax = self.canvas.axes.get_ylim()
        self.viewAllCoord = [xmin, xmax, ymin, ymax]
    else:
      xmin, xmax = self.canvas.axes.get_xlim()
      ymin, ymax = self.canvas.axes.get_ylim()
      self.viewAllCoord = [xmin, xmax, ymin, ymax]

  def viewAll(self):
    """
    (None) -> None
    Set new boundary values to the canvas according to values from the list 'viewAllCoord'.
    """
    try:
      gapX = abs(self.viewAllCoord[0] - self.viewAllCoord[1]) * 0.015
      gapY = abs(self.viewAllCoord[2] - self.viewAllCoord[3]) * 0.015
      self.canvas.axes.set_xlim(self.viewAllCoord[0] - gapX,
                                self.viewAllCoord[1] + gapX)
      self.canvas.axes.set_ylim(self.viewAllCoord[2] - gapY,
                                self.viewAllCoord[3] + gapY)
      self.canvas.draw()
      self.status.showMessage(
          "REMINDER: the 'Save' button records data from the graph according 'Axes range option' values, the 'View all' button does not change them.",
          3456)
    except (IndexError) as e:
      self.status.showMessage(
          'ERROR: failed to resize. Please check if there is some data plotted.',
          4567)

  def setCanvasBoundaries(self, x, X, y, Y):
    """
    (Float, Float, Float, Float) -> None
    Set new values to the fields from 'Axes range options'.
    """
    try:
      xmin = str(round(x, 3))
      xmax = str(round(X, 3))
      ymin = str(round(y, 3))
      ymax = str(round(Y, 3))

      self.plotXIntervalMin.setText(xmin)
      self.plotXIntervalMax.setText(xmax)
      self.plotYIntervalMin.setText(ymin)
      self.plotYIntervalMax.setText(ymax)
    except (IndexError, ValueError):
      self.status.showMessage('PASS.', 3456)

  def applyCanvasBoundaries(self):
    """
    (None) -> None
    Set new boundary values to the canvas according to the fields from 'Axes range options'.
    """
    try:
      pltYMin = float(self.plotYIntervalMin.text())
      pltYMax = float(self.plotYIntervalMax.text())
      pltXMin = float(self.plotXIntervalMin.text())
      pltXMax = float(self.plotXIntervalMax.text())

      if pltYMin < pltYMax:
        self.canvas.axes.set_ylim(pltYMin, pltYMax)
      else:
        self.status.showMessage(
            'Minimum Y-axis > Maximum Y-axis. Default Y-axis range applied.',
            3000)
      if pltXMin < pltXMax:
        self.canvas.axes.set_xlim(pltXMin, pltXMax)
      else:
        self.status.showMessage(
            'Minimum X-axis > Maximum X-axis. Default X-axis range applied.',
            3000)

      self.canvasInfo['user boundaries'] = [pltXMin, pltXMax, pltYMin, pltYMax]
      self.canvas.draw()
    except (IndexError, ValueError) as e:
      self.status.showMessage('PASS.', 3456)

  def changeSimulationOutput(self):
    """
    (None) -> None
    When eij files are loaded, update the data menus according to the selected simulation.
    """
    menus = [self.dataXMenu, self.dataYMenu]
    i = int(self.eijMenu.currentIndex())
    self.EijLabel = i + 1
    for i in menus:
      i.setDisabled(True)
      i.clear()
      i.setDisabled(False)
      i.addItems(self.labels)

  def applyMolDim(self):
    """
    (None) -> None
    Change the value of the RDF function according to the molecular dimensions entered by the user and show the new RDF function.
    """
    r = self.dataSet.loc[:, (self.GrLabel, 'r')]
    gr = self.dataSet.loc[:, (self.GrLabel, 'G(r)')]
    newGr = np.zeros(len(gr))

    gid = self.grMenu.currentIndex()
    graphHold = self.checkOverplot.isChecked()
    a = float(self.molDimA.text())
    b = float(self.molDimB.text())
    c = float(self.molDimC.text())

    halfDr = float(r[1] - r[0])
    esfConst = ((4 * 3.1415) / 3)

    for i, v in enumerate(gr):
      rMinus = r[i] - halfDr
      rPlus = r[i] + halfDr
      vEsf = (esfConst) * (rPlus**3 - rMinus**3)
      vPara = ((a + 2 * rPlus) * (b + 2 * rPlus) *
               (c + 2 * rPlus)) - ((a + 2 * rMinus) * (b + 2 * rMinus) *
                                   (c + 2 * rMinus))
      newGr[i] = ((2 * v * vEsf) / (vPara))

    self.setCanvasBoundaries(r.min(), r.max(), newGr.min(), newGr.max())

    self.integralMarker = False
    self.canvas.axes.axis('off')
    self.canvas.axes.set_position([0.075, 0.075, 0.7, 0.88])

    if not graphHold:
      self.canvas.axes.clear()
    self.canvas.axes.plot(r, newGr, linestyle='-')
    self.molDimAdjust = True
    if (self.horizontalGrLine == False or graphHold == False):
      self.canvas.axes.axhline(y=1, c='k', linestyle='-', label='_nolegend_')
      self.horizontalGrLine = True

    grLabel = self.grMenuItems[gid]
    self.plotTitles.append(r'${}*$'.format(grLabel))
    self.canvasInfo['title'] = '{}'.format(grLabel)
    self.canvas.axes.legend(self.plotTitles,
                            bbox_to_anchor=(1.05, 1),
                            loc=2,
                            ncol=1,
                            mode='expand',
                            borderaxespad=0.,
                            frameon=False,
                            numpoints=1,
                            prop={'size': 10},
                            handlelength=0.6,
                            borderpad=-0.8)

    self.canvasInfo['type'] = 'gr'
    self.canvasInfo['user data'] = newGr
    self.canvasInfo['user parameters'] = [a, b, c]

    self.updatePlotLimits()
    self.canvas.axes.axis('on')
    self.canvas.draw()

  def clickIntegral(self, event):
    """
    (MouseEvent) -> None
    Show an x-shaped symbol on the RDF curve and show the integral value corresponding to the clicked region.
    """
    if (self.extension == 'gr' and self.canvasInfo['type'] != 'histogram'):
      nPlots = len(self.canvas.axes.lines)
      if (self.checkOverplot.isChecked() == False):
        try:
          index = int(
              np.round((event.xdata - self.dataSet.loc[0,
                                                       (self.GrLabel, 'r')]) /
                       (self.dataSet.loc[1, (self.GrLabel, 'r')] -
                        self.dataSet.loc[0, (self.GrLabel, 'r')])))
          integral = self.nrData[index]
          self.integralLabel.setText(
              '<small><b>INTEGRAL OF G(R):</b></small>   N(%.3f) = %.3f' %
              (self.dataSet.loc[index, (self.GrLabel, 'r')], integral))

          if (self.integralMarker == False):
            self.canvas.axes.plot(event.xdata,
                                  event.ydata,
                                  'rx',
                                  label='_nolegend_')
            self.integralIndex = nPlots
            self.integralMarker = True
          else:
            self.canvas.axes.lines.pop(self.integralIndex)
            self.canvas.axes.plot(event.xdata,
                                  event.ydata,
                                  'rx',
                                  label='_nolegend_')

          self.canvas.draw()
        except (TypeError, IndexError) as e:
          self.status.showMessage('Please, click inside the graph.', 3000)
      else:
        self.status.showMessage(
            'Please, disable the OVERPLOT checkbutton to view the integral of the last plotted function.',
            3456)

  def savePlotData(self):
    """
    (None) -> None
    Write the values and information from the dictionary canvasInfo to a dat file. 
    """
    try:
      root = os.getcwd()
      timetag = time.strftime('%Y%m%d%H%M%S')
      fname = self.filename.split('/')[-1]
      title = ((self.canvasInfo['title']).replace('/', '')).replace('_', '')

      xmin, xmax, ymin, ymax = [
          float(i) for i in self.canvasInfo['user boundaries']
      ]
      #epsx = abs(xmax - xmin)*0.05
      #epsy = abs(ymax - ymin)*0.05
      #xmin, xmax, ymin, ymax = [xmin - epsx, xmax + epsx, ymin - epsy, ymax + epsy]

      # Line and Scatter
      if (self.canvasInfo['type'] == 'line') or (self.canvasInfo['type']
                                                 == 'scatter'):
        name = '{}_{}_{}_{}_{}.dat'.format(fname, self.extension,
                                           self.canvasInfo['type'], title,
                                           timetag)
        filename = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save file', name)[0]
        f = open('{}'.format(str(filename)), 'w')

        x = self.canvasInfo['data'].loc[:, self.Xlabel]
        y = self.canvasInfo['data'].loc[:, self.Ylabel]
        xmin, xmax, ymin, ymax = [
            float(i) for i in self.canvasInfo['user boundaries']
        ]

        if (len(self.eijMenuItems) != 0):
          f.write('# {}\n'.format(self.eijMenu.currentText()))

        title = (self.canvasInfo['title']).split('-')
        f.write('# {:>15}{:>15}\n'.format(title[0], title[1]))

        if (self.extension == 'gr'):
          Nr = self.nrData
          f.write('# {:>15}{:>20}{:>20}\n'.format('r',
                                                  self.dataYMenu.currentText(),
                                                  'N(r)'))
          for index, value in enumerate(x):
            if (value >= xmin) and (value <= xmax):
              if (y[index] >= ymin) and (y[index] <= ymax):
                f.write('{:>20e}{:>20e}{:20e}\n'.format(value, y[index],
                                                        Nr[index]))
        else:
          for index, value in enumerate(x):
            if (value >= xmin) and (value <= xmax):
              if (y[index] >= ymin) and (y[index] <= ymax):
                f.write('{:>20e}{:>20e}\n'.format(value, y[index]))

      # Histogram
      elif (self.canvasInfo['type'] == 'histogram'):
        name = '{}_{}_{}_{}_{}.dat'.format(fname, self.extension,
                                           self.canvasInfo['type'], title,
                                           timetag)
        filename = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save file', name)[0]
        f = open('{}'.format(str(filename)), 'w')

        bins = self.canvasInfo['data'][1]
        frequency = self.canvasInfo['data'][0]
        gaussian = self.canvasInfo['user data']
        xmin, xmax, ymin, ymax = [
            float(i) for i in self.canvasInfo['user boundaries']
        ]

        if (len(self.eijMenuItems) != 0):
          f.write('# {}\n'.format(self.eijMenu.currentText()))

        f.write('# {}\n'.format(self.canvasInfo['title']))
        f.write('# MEAN:      {:>15e}\n# STD. DEV.: {:>15e}\n'.format(
            self.canvasInfo['user parameters'][0],
            self.canvasInfo['user parameters'][1]))
        f.write('# {:>15}{:>20}{:>20}\n'.format('x', 'frequency(x)',
                                                'gaussian(x)'))
        for index, value in enumerate(bins):
          if (value >= xmin) and (value <= xmax):
            f.write('{:>20e}{:>20e}{:>20e}\n'.format(value, frequency[index],
                                                     gaussian[index]))

      # Autocorrelation
      elif (self.canvasInfo['type'] == 'autocorrelation'):
        name = '{}_{}_{}_{}_{}.dat'.format(fname, self.extension,
                                           self.canvasInfo['type'], title,
                                           timetag)
        filename = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save file', name)[0]
        f = open('{}'.format(str(filename)), 'w')

        if (len(self.eijMenuItems) != 0):
          f.write('# {}\n'.format(self.eijMenu.currentText()))

        f.write('# {} ({})\n'.format(self.canvasInfo['user parameters'][2],
                                     self.canvasInfo['title']))
        f.write('# USER GUESS\n')
        for index, value in enumerate(['A1', 'B1', 'A2', 'B2', 'A3', 'B3']):
          f.write('# {} = {}\n'.format(
              value, self.canvasInfo['user parameters'][0][index]))
        f.write('# BEST FIT\n')
        for index, value in enumerate(['A1', 'B1', 'A2', 'B2', 'A3', 'B3']):
          if index >= len(self.canvasInfo['user parameters'][1]):
            f.write('# {} = None\n'.format(value))
          else:
            f.write('# {} = {}\n'.format(
                value, self.canvasInfo['user parameters'][1][index]))
        f.write('# LARGEST T: {}\n'.format(
            self.canvasInfo['user parameters'][3]))

        t = self.canvasInfo['user data'][0]
        ct = self.canvasInfo['user data'][1]
        initial = self.canvasInfo['user data'][2]
        bestfit = self.canvasInfo['user data'][3]

        f.write('# {:>15}{:>20}{:>20}{:>20}\n'.format('t', 'C(t)', 'Initial',
                                                      'Best fit'))
        for index, value in enumerate(t):
          if (value >= xmin) and (value <= xmax):
            if (ct[index] >= ymin) and (ct[index] <= ymax):
              f.write('{:>20e}{:>20e}{:>20e}{:>20e}\n'.format(
                  value, ct[index], initial[index], bestfit[index]))

      # rdf
      elif (self.canvasInfo['type'] == 'gr'):
        name = '{}_{}_{}_{}_{}.dat'.format(fname, self.extension,
                                           self.canvasInfo['type'], title,
                                           timetag)
        filename = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save file', name)[0]
        f = open('{}'.format(str(filename)), 'w')

        f.write('# MOLECULAR DIMENSIONS\n')
        for index, value in enumerate(['A', 'B', 'C']):
          f.write('# {} = {}\n'.format(
              value, self.canvasInfo['user parameters'][index]))

        r = self.canvasInfo['data'][0]
        Gr = self.canvasInfo['data'][1]
        Nr = self.canvasInfo['data'][2]
        NewGr = self.canvasInfo['user data']

        ##print(xmin, xmax, ymin, ymax)

        f.write('# {:>15}{:>20}{:>20}{:>20}\n'.format('r', 'New G(r)', 'G(r)',
                                                      'N(r)'))
        for index, value in enumerate(r):
          if (value >= xmin) and (value <= xmax):
            if (NewGr[index] >= ymin) and (NewGr[index] <= ymax):
              f.write('{:>20e}{:>20e}{:>20e}{:>20e}\n'.format(
                  value, NewGr[index], Gr[index], Nr[index]))

      f.close()
      self.status.showMessage('File successfully saved.', 3456)

    except (IndexError, ValueError) as e:
      #print(e)
      self.status.showMessage(
          'ERROR: failed to save data. Please try to adjust the data intervals.',
          3456)
    except OSError as e:
      #print(e)
      self.status.showMessage('ERROR: failed to save data.', 3456)


def repeated_label(labels):
  lab_set = frozenset(labels)
  rep_ind = {lab: [] for lab in lab_set}

  for slab in lab_set:
    for i, lab in enumerate(labels):
      if lab == slab:
        rep_ind[slab].append(i)

  return rep_ind


def gen_new_labels(labels):
  reap = repeated_label(labels).values()
  for rinds in reap:
    j = 1
    for i in rinds:
      if len(rinds) > 1:
        labels[i] += str(j)
        j += 1

  return labels


def main():
  """
  (None) -> None
  Initialize and displace the window.
  """
  application = QtWidgets.QApplication(sys.argv)
  mainWindow = graphMainWindow()

  mainWindow.move(150, 50)

  mainWindow.show()
  sys.exit(application.exec_())


if __name__ == '__main__':
  main()
