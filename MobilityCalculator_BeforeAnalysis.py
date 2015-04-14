# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 10:05:31 2015

@author: Martin Friedl
"""

# for command-line arguments
import sys
import datetime
import re
# Python Qt4 bindings for GUI objects
from PyQt4 import QtCore,QtGui, uic
# Numpy functions for image creation
import numpy as np
# Pandas for manipulating datasets
import pandas as pd
# Matplotlib Figure object
#from matplotlib.figure import Figure
# import the Qt4Agg FigureCanvas object, that binds Figure to
# Qt4Agg backend. It also inherits from QWidget
#from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

class MyWindow(QtGui.QMainWindow):
    def __init__(self):
        super(MyWindow, self).__init__()
        uic.loadUi('MobilityCalculator.ui', self)
        
        self.initVariables()
        
        # Connect up the buttons.
        self.pbLoad.clicked.connect(self.loadDataClickedSlot)
        self.pbSavePltDat.clicked.connect(self.savePltDat)   
        self.cbXAxis.currentIndexChanged.connect(self.xAxisChangedSlot)
        self.cbYAxis.currentIndexChanged.connect(self.yAxisChangedSlot)
        self.cbLineStyle.currentIndexChanged.connect(self.lineStyleChangedSlot)      
        self.cbMarkerStyle.currentIndexChanged.connect(self.markerStyleChangedSlot)
        self.cbPlotColour.currentIndexChanged.connect(self.plotColourChangedSlot)        
        self.ckbAvgData.stateChanged.connect(self.avgDataChangedSlot)
        self.ckbRmvFiller.stateChanged.connect(self.rmvFillerChangedSlot)  
        self.tabWidget.currentChanged.connect(self.tabChangedSlot)
        
        #Add NavigationToolbar
        self.mplwidget.figure.set_dpi(120)
        self.mplwidget.mpl_toolbar = NavigationToolbar(self.mplwidget.figure.canvas, self.mplwidget)
        self.verticalLayout.setDirection(QtGui.QBoxLayout.BottomToTop)
        self.verticalLayout.addWidget(self.mplwidget.mpl_toolbar,1)

        self.mplwidget.myplot = self.mplwidget.figure.add_subplot(111)        
        self.mplwidget.myplot.set_title("Pulse Curve")
        self.mplwidget.myplot.set_xlabel("Pt Index")
        self.mplwidget.myplot.set_ylabel("Voltage (V)")  
        
        self.show()
        self.updatePlot()    
        
        self.lineStyleDict = {0:'-',1:'--',2:'-.',3:':',4:''}
        self.markerDict = {0:'.',1:'o',2:'+',3:'.',4:'1',5:'2',6:'3',7:'4',8:''}
        self.lineColourDict = {0:'b',1:'g',2:'r',3:'c',4:'m',5:'y',6:'k',7:'w'}

        self.dict_NWLen = {'MF-04-IA' : {'AD':3.400E-6, 'AB':0.924E-6, 'BC':0.924E-6,'CD':0.924E-6},
                           'MF-03-IA' : {'AD':2.971E-6, 'AB':None, 'BC':None,'CD':None}}


    def initVariables(self):
        
        #Filename Information        
        self.wafer = None
        self.NWType = None
        self.deviceID = None
        self.measurement = None
        self.lighting = None
        self.testType = None
#        self.infoDict = None #Dictionary linking table item names to row indexes

        #Header File Information
        self.delayInit = None
        self.delayAfterSet = None
        self.delayAfterRead = None                 
        self.numPoints = None                 
        self.numMeasurements = None
        self.numSetpoints = None

        #DataSet Initialization
        self.headers = None
        self.units = None
        self.dataArray_raw = None
        self.dataArray_avg = None
        self.dataArray_filt = None        
        self.dataArray_avgFilt = None
        self.pltXData = None
        self.pltYData = None

        self.cbLineStyle.currentIndexChanged.connect(self.lineStyleChangedSlot)      
        self.cbMarkerStyle.currentIndexChanged.connect(self.markerStyleChangedSlot)
        self.cbPlotColour.currentIndexChanged.connect(self.plotColourChangedSlot)

    def tabChangedSlot(self):
        selected_index = self.tabWidget.currentIndex
        if selected_index == 0: #Plot chooser
            pass
        elif selected_index == 1: #Analysis
            pass
        elif selected_index == 2: #Calculated Values
            pass       

    def avgDataChangedSlot(self):               
        self.updatePlot()
    def rmvFillerChangedSlot(self):
        self.updatePlot()
    def lineStyleChangedSlot(self):
        self.mplwidget.myplot.lines[0].set_linestyle(self.lineStyleDict[self.cbLineStyle.currentIndex()])
        self.mplwidget.figure.canvas.draw()          
    def markerStyleChangedSlot(self):
        self.mplwidget.myplot.lines[0].set_marker(self.markerDict[self.cbMarkerStyle.currentIndex()])
        self.mplwidget.figure.canvas.draw() 
    def plotColourChangedSlot(self):
        self.mplwidget.myplot.lines[0].set_color(self.lineColourDict[self.cbPlotColour.currentIndex()])
        self.mplwidget.figure.canvas.draw() 
    def xAxisChangedSlot(self):
        self.changeXAxisData()
    def yAxisChangedSlot(self):
        self.changeYAxisData()        
    def changeXAxisData(self):
        self.updatePlot()
    def changeYAxisData(self):
        self.updatePlot()        
    def loadDataClickedSlot(self): #Wrapper function since must have compatibility between signal and slot
        self.loadData()
    def loadData(self,fileName = None):
        if fileName is None:
            fileName = QtGui.QFileDialog.getOpenFileName(self, 'Choose data file to open', '.', filter='*.lvm')
        if fileName:
            if self.dataArray_raw is not None: 
                self.initVariables()
            self.filename = fileName
            self.parseFilename()
            f = open(self.filename, 'r')
            self.readFileHeader(f)
            f.seek(0) #Reset cursor to beginning of the file
            self.dataArray_raw = self.importData(f)
            f.close()          
            self.populateComboBoxes() #Populate plot choosing combo boxes         

            #Create various filtered/averaged data arrays for viewing
            self.dataArray_avg=self.averageData(self.dataArray_raw) #Average data over number of data points at each setpoint
            self.dataArray_filt=self.removeFillerPoints(self.dataArray_raw,duty=8,ppS=self.numPoints)
            self.dataArray_avgFilt=self.removeFillerPoints(self.dataArray_avg,duty=8,ppS=1) 
                
            #Put data into a pandas dataframe for later use
            self.dataFrame = pd.DataFrame(data=zip(self.dataArray_avgFilt[:,0],self.dataArray_avgFilt[:,1],self.dataArray_avgFilt[:,2],self.dataArray_avgFilt[:,3]),columns=[self.headers[0],self.headers[1],self.headers[2],self.headers[3]])          
            self.dataFrame = pd.pivot_table(self.dataFrame,values='Is', index='Vds', columns='Vgs')
            print "Loaded file: " + fileName
        self.updatePlot()
    
    #Averages the data over the number of data points taken for each setpoint
    def averageData(self,dataSet):
        columns = np.shape(dataSet)[-1]
        rows = np.shape(dataSet)[0]
        dataSet = np.reshape(dataSet,(rows/self.numPoints,self.numPoints,columns))        
        dataSet = dataSet.mean(1)
        return dataSet
    
    #Removes all of the points in between the data points that we want to analyze
    def removeFillerPoints(self,dataSet,duty=8,ppS=1):
        columns = np.shape(dataSet)[-1]
        rows = np.shape(dataSet)[0]
        period = duty+2 #One extra point for the data point and the second for the negative pulse
        dataSet = np.reshape(dataSet,(rows/period/ppS,period*ppS*columns))
        dataSet = dataSet[:,duty*ppS*columns:(duty+1)*ppS*columns]   
        dataSet = np.reshape(dataSet,(rows/(duty+2),columns))
        return dataSet
    
    def importData(self,f):
        self.headers = None
        self.units = None
        dataArray = None
        lineDataList = []
        for line in f: 
            if (line[0] == '#')or(line == '\n'): continue
            splitLine = line[:-1].split('\t')
            if self.headers is None:
                self.headers = splitLine
                continue
            if self.units is None:
                self.units = splitLine
                continue  
            try: #Try to convert line to list of floats
                lineData = map(float,splitLine)
            except: #If it fails, we have strings on that line so skip to next line
                continue
            if not(len(lineDataList)==0): #If the array is not empty
                if len(lineDataList[0])!=len(map(float,splitLine)): #Make sure the line and the array have the same number of columns
                    continue                    
            lineDataList.append(lineData)       
        dataArray=np.vstack(lineDataList) #Stack list up into an array
        dataArray[dataArray > 999999999] = np.nan #Remove overflow values               
        return dataArray
            
    def populateComboBoxes(self):
        self.cbXAxis.clear()
        self.cbXAxis.addItems(self.headers)
        self.cbYAxis.clear()        
        self.cbYAxis.addItems(self.headers)
        
    def readFileHeader(self,f):      
        for line in f:
            if line == '\n': #Skip empty lines
                continue
            if line[0] != "#": #Stop if we reach the end of the header
                break
            if "# Delay(init):" in line:
                parsedLine = line.split(',')
                for text in parsedLine:
                    splitText = text.split(":")
                    value = splitText[-1].strip() #Strip whitespace
                    value = int(re.sub(r'[^0-9]+','',value)) #Remove non-integer characters
                    if "Delay(init)" in text:
                        self.delayInit = value
                    elif "delay(after set)" in text:
                        self.delayAfterSet = value
                    elif "delay(bef.read)" in text:
                        self.delayAfterRead = value                    
                    elif "np. of points" in text:
                        self.numPoints = value                   
                    elif "np. of measurements" in text:
                        self.numMeasurements = value                   
                    elif "np. of setpoints" in text:
                        self.numSetpoints = value             
                    else:
                        pass
        self.updateHeaderInfoTable()

    def parseFilename(self): #Extract information from the data filename
        parsed = str(self.filename.split("/")[-1]).split(".")[0].split("_")
        self.wafer = parsed[0]
        self.NWType = parsed[1]
        self.deviceID = parsed[2].split("-")[0]+"-"+parsed[2].split("-")[1]
        self.measurement = parsed[2].split("-")[-1]
        self.lighting = parsed[3]
        self.testType = parsed[4]
        self.updateDeviceInfoTable()
        
    def updateDeviceInfoTable(self): #Updates table with information from filename of data file
        self.tblDeviceInfo.item(self.findItemIndexes(self.tblDeviceInfo,"Wafer ID")[0],1).setText(self.wafer)
        self.tblDeviceInfo.item(self.findItemIndexes(self.tblDeviceInfo,"Nanowire ID")[0],1).setText(self.NWType)
        self.tblDeviceInfo.item(self.findItemIndexes(self.tblDeviceInfo,"Device ID")[0],1).setText(self.deviceID)
        self.tblDeviceInfo.item(self.findItemIndexes(self.tblDeviceInfo,"Measurement #")[0],1).setText(self.measurement)
        self.tblDeviceInfo.item(self.findItemIndexes(self.tblDeviceInfo,"Lighting")[0],1).setText(self.lighting)
        self.tblDeviceInfo.item(self.findItemIndexes(self.tblDeviceInfo,"Test Type")[0],1).setText(self.testType) 
    
    def updateHeaderInfoTable(self): #Updates table with information from header of data file
        self.tblHeaderInfo.item(self.findItemIndexes(self.tblHeaderInfo,"Initial Delay")[0],1).setText(str(self.delayInit)+"ms")
        self.tblHeaderInfo.item(self.findItemIndexes(self.tblHeaderInfo,"Delay After Set")[0],1).setText(str(self.delayAfterSet)+"ms")
        self.tblHeaderInfo.item(self.findItemIndexes(self.tblHeaderInfo,"Delay After Read")[0],1).setText(str(self.delayAfterRead)+"ms")
        self.tblHeaderInfo.item(self.findItemIndexes(self.tblHeaderInfo,"Num Points")[0],1).setText(str(self.numPoints))
        self.tblHeaderInfo.item(self.findItemIndexes(self.tblHeaderInfo,"Num Measurements")[0],1).setText(str(self.numMeasurements))
        self.tblHeaderInfo.item(self.findItemIndexes(self.tblHeaderInfo,"Num Setpoints")[0],1).setText(str(self.numSetpoints))          
    
    #Find index of row header in a QTable
    def getRowHeaderIndex(self,table,searchString):
        output=None        
        for i in range(table.rowCount()):
            if self.tblDeviceInfo.verticalHeaderItem(i).text() == searchString:
                output = i
                break
        return output
    #Find index of column header in a QTable
    def getColumnHeaderIndex(self,table,searchString):
        output=None        
        for i in range(table.rowCount()):
            if self.tblDeviceInfo.horizontalHeaderItem(i).text() == searchString:
                output = i
                break
        return output
    #Find the x,y indices of an item in a QTable
    def findItemIndexes(self,table,searchString):
        output=np.array([])
        filteredItems = table.findItems(searchString,QtCore.Qt.MatchExactly)        
        for item in filteredItems:
            output = np.append(output,(item.row(),item.column()))
        return output                

    def savePltDat(self):
        fileName = QtGui.QFileDialog.getSaveFileName(self, 'Choose save data file name', '.', filter='*.dat')
        f = open(fileName, 'w')
        f.write("# Created by MobilityCalculator - Martin Friedl\n")
        f.write("# "+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+"\n")
        f.write('# '+self.headers[self.cbXAxis.currentIndex()]+" ("+self.units[self.cbXAxis.currentIndex()]+")"+'\t'+self.headers[self.cbYAxis.currentIndex()]+" ("+self.units[self.cbYAxis.currentIndex()]+")"+'\n')
        for x,y in zip(self.pltXData,self.pltYData):
            f.write('%E\t' % x +'%E\t' % y + '\n')
        f.close()            
        print "Saved to " + str(self.filename)
        
    def updatePlot(self):
        
        try:
            if self.ckbAvgData.checkState() and not(self.ckbRmvFiller.checkState()):
                self.pltXData = self.dataArray_avg[:,self.cbXAxis.currentIndex()]
                self.pltYData = self.dataArray_avg[:,self.cbYAxis.currentIndex()]            
            elif not(self.ckbAvgData.checkState()) and self.ckbRmvFiller.checkState():
                self.pltXData = self.dataArray_filt[:,self.cbXAxis.currentIndex()]
                self.pltYData = self.dataArray_filt[:,self.cbYAxis.currentIndex()] 
            elif self.ckbAvgData.checkState() and self.ckbRmvFiller.checkState():
                self.pltXData = self.dataArray_avgFilt[:,self.cbXAxis.currentIndex()]
                self.pltYData = self.dataArray_avgFilt[:,self.cbYAxis.currentIndex()] 
            else:
                self.pltXData = self.dataArray_raw[:,self.cbXAxis.currentIndex()]
                self.pltYData = self.dataArray_raw[:,self.cbYAxis.currentIndex()]       
        except:
            return
        
        
#        if (self.pltXData is None) or (self.pltYData is None):
#            return
            
        self.mplwidget.myplot.cla()
        marker = self.markerDict[self.cbMarkerStyle.currentIndex()]
        linestyle = self.lineStyleDict[self.cbLineStyle.currentIndex()]
        colour = self.lineColourDict[self.cbPlotColour.currentIndex()]
        self.mplwidget.myplot.plot(self.pltXData,self.pltYData,marker=marker, linestyle=linestyle, color=colour)
        self.mplwidget.myplot.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        self.mplwidget.myplot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        self.mplwidget.myplot.set_xlabel(self.headers[self.cbXAxis.currentIndex()]+" ("+self.units[self.cbXAxis.currentIndex()]+")")
        self.mplwidget.myplot.set_ylabel(self.headers[self.cbYAxis.currentIndex()]+" ("+self.units[self.cbYAxis.currentIndex()]+")")
        self.mplwidget.myplot.set_title(self.wafer+" "+self.deviceID+" Exp:"+self.measurement)
        self.mplwidget.myplot.grid()
        self.mplwidget.figure.canvas.draw()          
        
if __name__ == '__main__':
    # Create the GUI application
    qApp = QtGui.QApplication(sys.argv)
    window = MyWindow()
    window.show() 
    window.loadData("D:/Dropbox/LMSC/Programming/ACMobility/MF-04-IA_04-23B_96-26-0_dark_MobilityPulsedAC-AD.lvm")
#    window.loadData("D:/Dropbox/LMSC/Programming/ACMobility/MF-04-IA_04-23B_92-26-1_dark_MobilityPulsedAC-AD.lvm")

    # start the Qt main loop execution, exiting from this script
    # with the same return code of Qt application
    sys.exit(qApp.exec_()) 
    
    