# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 10:05:31 2015

@author: Martin Friedl
"""

#TODO: Get the color map working
#TODO: Change to raw data and analysis data to a two-tab environment with separate MPL canvases
#TODO: Add effective mobility at 0V to the analysis results table

# for command-line arguments
import sys
import datetime
# Python Qt4 bindings for GUI objects
from PyQt4 import QtCore,QtGui, uic
# Numpy functions for image creation
import numpy as np
# Pandas for manipulating datasets
import pandas as pd
# Scipy Signal library for Savitzky-Golay filtering of data
import scipy.signal
# Matplotlib Figure object
#from matplotlib.figure import Figure
# import the Qt4Agg FigureCanvas object, that binds Figure to
# Qt4Agg backend. It also inherits from QWidget
#from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

from DataLoader import dataFile

class MyWindow(QtGui.QMainWindow):
    def __init__(self):
        super(MyWindow, self).__init__()
        uic.loadUi('MobilityCalculator.ui', self)
        
        self.initVariables()
        
        # Connect up the buttons.
        self.pbLoad.clicked.connect(self.loadDataClickedSlot)
        self.pbSavePltDat.clicked.connect(self.savePltDat)   
        self.cbXAxis.currentIndexChanged.connect(self.updatePlot_RawData)
        self.cbYAxis.currentIndexChanged.connect(self.updatePlot_RawData)        
        self.cbLineStyle.currentIndexChanged.connect(self.lineOrMarkerStyleChangedSlot)      
        self.cbMarkerStyle.currentIndexChanged.connect(self.lineOrMarkerStyleChangedSlot)
        self.cbPlotColour.currentIndexChanged.connect(self.plotColourChangedSlot)       
        self.ckbAvgData.stateChanged.connect(self.updatePlot_RawData)
        self.ckbRmvFiller.stateChanged.connect(self.updatePlot_RawData)  
        self.tabWidget.currentChanged.connect(self.tabChangedSlot)        
        self.cbAXAxis.currentIndexChanged.connect(self.updatePlotFromAxisSettings)
        self.cbAYAxis.currentIndexChanged.connect(self.updatePlotFromAxisSettings)        
        self.cbAIndex.currentIndexChanged.connect(self.updatePlotFromAxisSettings)
        self.sbLegendEntries.valueChanged.connect(self.updatePlotFromAxisSettings)
        self.cbLineStyle_Analysis.currentIndexChanged.connect(self.lineOrMarkerStyleChangedSlot)      
        self.cbMarkerStyle_Analysis.currentIndexChanged.connect(self.lineOrMarkerStyleChangedSlot)
        self.cbPlotColour_Analysis.currentIndexChanged.connect(self.AplotColourmapChangedSlot)  
        self.lwMeasSelect.itemSelectionChanged.connect(self.updateDataSet)
        self.lwSweepDirSelect.itemSelectionChanged.connect(self.updateDataSet)        
        self.ckbShowVtFit.stateChanged.connect(self.updatePlotFromAxisSettings)
        self.ckbSmoothing.stateChanged.connect(self.updateDataSet)
        self.hsSmoothing.valueChanged.connect(self.updateDataSet)
        self.ckbDriftCorr.stateChanged.connect(self.updateDataSet)
        
        self.leFitVtFrom.setValidator(QtGui.QDoubleValidator())
        self.leFitVtTo.setValidator(QtGui.QDoubleValidator())
        self.leFitVtFrom.editingFinished.connect(self.updateDataSet)         
        self.leFitVtTo.editingFinished.connect(self.updateDataSet)
        self.leFitVtFrom.setText(str(self.defVtFitFrom))
        self.leFitVtTo.setText(str(self.defVtFitTo))
        
        #Add NavigationToolbar
        self.mplwidget.figure.set_dpi(120)
        self.mplwidget.mpl_toolbar = NavigationToolbar(self.mplwidget.figure.canvas, self.mplwidget)
        self.verticalLayout.setDirection(QtGui.QBoxLayout.BottomToTop)
        self.verticalLayout.addWidget(self.mplwidget.mpl_toolbar,1)

        self.mplwidget.myplot = self.mplwidget.figure.add_subplot(111)        
        self.mplwidget.myplot.set_title("Pulse Curve")
        self.mplwidget.myplot.set_xlabel("Pt Index")
        self.mplwidget.myplot.set_ylabel("Voltage (V)")  
        
        self.pltCfg_RawData = {'marker':'.','linestyle':'-','linecolor':'b','title':'Title','xlabel':'xlabel','ylabel':'ylabel'}        
        self.pltCfg_Analysis = {'marker':'.','linestyle':'-','linecolor':'b','title':'Title','xlabel':'xlabel','ylabel':'ylabel'}                
        
        self.lineStyleDict = {0:'-',1:'--',2:'-.',3:':',4:''}
        self.markerDict = {0:'.',1:'o',2:'+',3:'.',4:'1',5:'2',6:'3',7:'4',8:''}
        self.lineColourDict = {0:'b',1:'g',2:'r',3:'c',4:'m',5:'y',6:'k',7:'w'}

        self.dict_NWLen = {'MF-08-ZA' : {'AD':2.700E-6, 'AB':0.929E-6, 'BC':0.929E-6,'CD':0.929E-6},
                           'MF-07-GA' : {'AD':3.400E-6, 'AB':0.924E-6, 'BC':0.924E-6,'CD':0.924E-6},
                           'MF-06-IA' : {'AD':2.667E-6, 'AB':0.745E-6, 'BC':0.745E-6,'CD':0.745E-6},
                           'MF-05-IA' : {'AD':3.400E-6, 'AB':0.924E-6, 'BC':0.924E-6,'CD':0.924E-6},
                           'MF-04-IA' : {'AD':3.400E-6, 'AB':0.924E-6, 'BC':0.924E-6,'CD':0.924E-6},                                                      
                           'MF-03-IA' : {'AD':2.971E-6, 'AB':None, 'BC':None,'CD':None}}
        self.dict_NWDiam = {'MF-08-ZA' : 83E-9,  #Zn3As2
                            'MF-07-GA' : 100E-9, #InAs/GaSb Core-shell
                            'MF-06-IA' : 130E-9, #InAsSb
                            'MF-05-IA' : 100E-9, #InAs 14-04-23B
                            'MF-04-IA' : 100E-9, #InAs 14-04-23B
                            'MF-03-IA' : 100E-9} #InAs 14-04-23B
        self.show()
        self.updatePlot_RawData() 

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
        self.dataArray_orig = None        
        self.dataArray_raw = None
        self.dataArray_avg = None
        self.dataArray_filt = None        
        self.dataArray_avgFilt = None
        self.dataArray_avgFilt_driftCorr = None        
        self.pltXData = None
        self.pltYData = None
        
        self.dataFrame_orig = None
        self.dataFrame = None        
        self.dataFrame_mobility = None
        self.dataFrame = None     
        
        #Nanowire parameters
        self.effMass = 0.023 #effective electron mass in gamma valley of InAs
        self.NWRadius = 50E-9
        self.NWLength = 3.400E-6
        self.oxideThickness = 200E-9
        self.corrDielcConst = 2.25 #adjusted for NW not being embedded in SiO2 [Wunnicke 2006]
      
        #Defaults
        self.defVtFitFrom = -5
        self.defVtFitTo = 0        
      
        #Analysis Results
        self.Vt = None
        self.VtErr = None
        self.peakMobility = None
        self.electronConcent = None   
        self.subthresholdSlope = None
        self.onOffRatio = None
        self.sweepRate = None
        self.rho = None
        self.mfp = None
        self.cohTime = None
        self.driftVelocity = None
        
        self.dataFile = None
      
    def calculateMobilty(self):
        q = 1.60217657E-19
        er = self.corrDielcConst 
        tox = self.oxideThickness
        radius = self.NWRadius
        nwLen = self.NWLength
        
        Cox = self.calculateOxideCapacitance(er,nwLen,radius,tox)
        
        Vt = []
        peakMobility = []
        electronConc = []
        subThresh = []
        rho = []
        onOffRatio = []
        sweepRate = []
        mfp = []
        cohTime = []
        driftV = []
        
        fitVtFrom = float(self.leFitVtFrom.text()) 
        fitVtTo = float(self.leFitVtTo.text()) 

        allDat = pd.DataFrame()        
        for Vds in self.dataFrame['Vds'].unique():
            if Vds == 0:
                continue
            dat = self.dataFrame.loc[self.dataFrame['Vds']==Vds].copy()
            #Average the duplicate data points for all possible Vgs and Vds pairs,
            dat = dat.reset_index(drop=True).groupby(["Vgs","Vds"]).mean().reset_index()        
            
#            dat['Is'].loc[(dat['Vgs']==-3) | (dat['Vgs']==-2) | (dat['Vgs']==-1)] = np.NaN #Erase faulty data points         
#            dat = dat.interpolate(method = 'cubic') #Linear interpolation of missing data points
            
            dat['Is_smooth'] = dat['Is']
            if (self.ckbSmoothing.isChecked()): #Smooth data with Savitzky-Golay filter
                smoothVal = self.hsSmoothing.value()*2+1
                dat['Is_smooth'] = scipy.signal.savgol_filter(dat['Is_smooth'],smoothVal,2) #Parameters are window size and polynomial order
            #Find the threshold voltage
            x = dat.loc[(dat['Vgs'] >= fitVtFrom) & (dat['Vgs'] <= fitVtTo)]['Vgs']
            y = dat.loc[(dat['Vgs'] >= fitVtFrom) & (dat['Vgs'] <= fitVtTo)]['Is_smooth']
            p = np.polyfit(x,y,1)
            dat['Is_Vtfit']= dat['Vgs']*p[0]+p[1]
            Vt.append(-p[1]/p[0])
    
            #Calculate the conductance
            dat['Transconductance'] = (dat['Is_smooth'].copy().diff().abs())/(dat['Vgs'].diff().abs())
            if (self.ckbSmoothing.isChecked()): #Smooth conductance curve with Savitzky-Golay filter
                smoothVal = self.hsSmoothing.value()*2+1
                dat['Transconductance'] = scipy.signal.savgol_filter(dat['Transconductance'],5,3) #Parameters are window size and polynomial order           

            #Calculate Field Effect Mobility
            dat['Mobility_fe'] = np.abs(np.power(nwLen,2)*dat['Transconductance']/Cox/Vds)
            peakMobility.append(dat['Mobility_fe'].max())
            #Calculate Charge Concentration - From Sourribes 2014
            dat['n_electrons'] = Cox*(dat['Vgs']-Vt[-1])/(q*np.pi*np.power(radius,2)*nwLen)
            electronConc.append(dat.loc[dat['Vgs'] == 0]['n_electrons'].mean())
            #Calculate subthreshold slope
            fitSubThreshFrom = dat['Vgs'].min()
            fitSubThreshTo = Vt[-1]
            try:                
                p = np.polyfit(dat.loc[(dat['Vgs'] > fitSubThreshFrom) & (dat['Vgs'] < fitSubThreshTo)]['Vgs'],np.log10(dat.loc[(dat['Vgs'] > fitSubThreshFrom) & (dat['Vgs'] < fitSubThreshTo)]['Is_smooth']),1)
                subThresh.append(1/p[0])
            except:
                subThresh = np.nan #Didn't work
            #Calculate On/Off Ratio
            onOffRatio.append(np.abs(dat['Is_smooth'].max()/dat['Is_smooth'].min()))
            #Calculate Gate Voltage Sweep Rate
            sweepRate = (dat['Vgs'].diff().abs()/dat['t'].diff().abs()).mean()
            #Calculate Nanowire Resistivity
            rho.append(self.calcResistivity(np.abs(Vds/dat.loc[dat['Vgs'] == 0]['Is_smooth']).mean(),radius,nwLen))
            #Calculate Electron Mean Free Path
            a, b, c = self.calcMeanFreePath(electronConc[-1],dat.loc[dat['Vgs'] == 0]['Is_smooth'].mean(),dat.loc[dat['Vgs'] == 0]['Mobility_fe'].mean(),radius,self.effMass)     
            mfp.append(a)
            cohTime.append(b)
            driftV.append(c)
            allDat=allDat.append(dat,ignore_index=True)
        
        self.Vt = np.mean(Vt)        
        self.VtErr = np.std(Vt,ddof=1)/np.sqrt(len(Vt)) #Standard error of sample mean
        if len(peakMobility)>2:
            peakMobility = self.reject_outliers(np.array(peakMobility))
        self.peakMobility = np.mean(peakMobility)
        self.electronConcent = np.mean(electronConc)
        
        self.subthresholdSlope = np.mean(subThresh)
        self.rho = np.mean(rho)
        self.onOffRatio = np.mean(onOffRatio)
        self.sweepRate = np.mean(sweepRate)
        self.mfp = np.mean(mfp)
        self.cohTime = np.mean(cohTime)
        self.driftVelocity = np.mean(driftV)
        
        allDat['ChanConductance'] = np.nan
        allDat['EffMobility'] = np.nan
        
        # Calculate effective mobility
        for Vgs in allDat['Vgs'].unique():
            dat = allDat.loc[allDat['Vgs']==Vgs].copy()
            x = dat['Vds']
            y = dat['Is_smooth']
            p = np.polyfit(x,y,1)
            allDat.loc[allDat['Vgs']==Vgs,'ChanConductance'] = p[0]
            allDat.loc[allDat['Vgs']==Vgs,'EffMobility'] = np.abs(np.power(nwLen,2)*p[0]/Cox/(Vgs-self.Vt))
        
        self.dataFrame=allDat       
        self.updateAnalysisInfoTables()
    
    def reject_outliers(self,data,m=2.):
        d = np.abs(data - np.median(data))
        mdev = np.median(d)
        s = d/mdev if mdev else 0.
        return data[s<m]    
    
    #Calculate the mean free path
    def calcMeanFreePath(self,n,I,mu,r,m_eff):
        A = np.pi*np.power(r,2) #Nanowire cross-sectional area
        e = 1.602E-19 #Elementary charge
        me = 9.109E-31 #Electron mass
        m = me*m_eff #Effective electron mass
        vd = I/A/e/n #Drift velocity
        t_coh = mu*m/e #Coherence time
        mfp = vd*t_coh; #Mean free path 
        return (mfp, t_coh,vd)
    
    #Get resistivity of a nanowire
    def calcResistivity(self,R,r,l):        
        csA = np.pi*np.power(r,2)
        return R*csA/l
        
    #Semi-empirical calculation of capacitance of NW on infinite gate plane
    def calculateOxideCapacitance(self,er,l,r,tox):
        e0 = 8.854187817E-12 
        C_ox = 2.0*np.pi*e0*er*l/np.arccosh((r+tox)/r);
        return C_ox        
        
    def tabChangedSlot(self):
        selected_index = self.tabWidget.currentIndex()
        if selected_index == 0: #Plot chooser
            self.updatePlot_RawData()
        elif selected_index == 1: #Analysis
            self.updatePlotFromAxisSettings()
    def lineOrMarkerStyleChangedSlot(self):
        for l in self.mplwidget.myplot.lines:
            if self.tabWidget.currentIndex() == 0:
                l.set_linestyle(self.lineStyleDict[self.cbLineStyle.currentIndex()])
                l.set_marker(self.markerDict[self.cbMarkerStyle.currentIndex()])
            else:
                l.set_linestyle(self.lineStyleDict[self.cbLineStyle_Analysis.currentIndex()])                
                l.set_marker(self.markerDict[self.cbMarkerStyle_Analysis.currentIndex()])
        self.mplwidget.figure.canvas.draw()        

    def AplotColourmapChangedSlot(self):
        print "Not working..."

    def plotColourChangedSlot(self):
        self.mplwidget.myplot.lines[0].set_color(self.lineColourDict[self.cbPlotColour.currentIndex()])
        self.mplwidget.figure.canvas.draw()        
        
    def loadDataClickedSlot(self): #Wrapper function since must have compatibility between signal and slot
        fileName = QtGui.QFileDialog.getOpenFileName(self, 'Choose data file to load', self.dataFile.fileInfo['Filename'], filter='*.dat;;*.lvm;;*.*')
        if fileName:
            self.loadData(fileName=fileName)

    def loadData(self,fileName = None):
        self.dataFile = dataFile(fileName)
        self.lwSweepDirSelect.item(0).setSelected(True) #Select forward sweep direction as default        
        self.setNWGeometry(self.dataFile.fileInfo['Wafer'],self.dataFile.fileInfo['Test Type'])
#        self.dataFile.duty = 8
        self.processRawData(self.dataFile.dataArray,self.dataFile.headerInfo['Num Points'],self.dataFile.headerInfo['Num Setpoints'])       
        self.calculateMobilty()   
        self.dataFile.updateDeviceInfoTable(self.tblDeviceInfo)
        self.dataFile.updateHeaderInfoTable(self.tblHeaderInfo)                
        self.populateComboBoxes(self.dataFile) #Populate plot choosing combo boxes         
        self.populateMeasurementList(self.dataFile) #Populate list to choose measurements to use            
        print "Loaded file: " + fileName
        self.updatePlot_RawData()
    
    def setNWGeometry(self,wafername,testType):
        contactsUsed = testType[-2:]
        try:
            self.NWLength =  self.dict_NWLen[wafername][contactsUsed]
            self.NWRadius = self.dict_NWDiam[wafername]/2.
        except:
            print 'ERROR: Unknown test type, using default NW geometry!'
            self.NWLength = 1E-6
            self.NWRadius = 100E-9
    
    #Create various filtered/averaged data arrays for viewing
    def processRawData(self,dataArray,ppSP,SPs):
        self.dataArray_avg=self.averageData(dataArray,ppSP) #Average data over number of data points at each setpoint
        self.dataArray_filt=self.removeFillerPoints(dataArray,ppSP=ppSP,duty=self.dataFile.duty,SPs=SPs)
        self.dataArray_avgFilt=self.removeFillerPoints(self.dataArray_avg,ppSP=1,duty=self.dataFile.duty,SPs=SPs) 
        self.dataArray_avgFilt_driftCorr = self.correctMeasDrift(self.dataArray_avg,self.dataArray_avgFilt,ppSP=1,duty=self.dataFile.duty,SPs=SPs)
        #For the analysis tab
        if self.ckbDriftCorr.isChecked():
            self.dataFrame_orig = pd.DataFrame(data=zip(self.dataArray_avgFilt_driftCorr[:,0],self.dataArray_avgFilt_driftCorr[:,1],self.dataArray_avgFilt_driftCorr[:,2],self.dataArray_avgFilt_driftCorr[:,3]),columns=[self.dataFile.headers[0],self.dataFile.headers[1],self.dataFile.headers[2],self.dataFile.headers[3]])              
        else:            
            self.dataFrame_orig = pd.DataFrame(data=zip(self.dataArray_avgFilt[:,0],self.dataArray_avgFilt[:,1],self.dataArray_avgFilt[:,2],self.dataArray_avgFilt[:,3]),columns=[self.dataFile.headers[0],self.dataFile.headers[1],self.dataFile.headers[2],self.dataFile.headers[3]])          
        self.dataFrame = self.dataFrame_orig
#        self.updateSmoothingSlider(5,self.oddNum(len(self.dataFrame['Vgs'].unique())),9)
        self.updateSmoothingSlider(5,15,5)

    #Over measurement time period, the current will drift. This function corrects for this drift looking at the 
    #relative CHANGE in voltage between the last Vgs=0V and each respective Vgs=x data point
    def correctMeasDrift(self,dat,datAvgFilt,ppSP=1,duty=8,SPs=1):
        return datAvgFilt
        period = duty+2
        tmp = dat[:,2] #Extract only the current column
        tmp = np.reshape(tmp,(-1,period))
        datAvgFilt[:,2] = tmp[:,period-2]-tmp[:,period-3] #Use only the differences from the baseline
        return datAvgFilt
            
    def oddNum(self,num):
        return (num+1 & ~1) - 1
    #Averages the data over the number of data points taken for each setpoint
    def averageData(self,dataArray,numPts):
        columns = np.shape(dataArray)[-1]
        rows = np.shape(dataArray)[0]
        dataSet = np.reshape(dataArray,(rows/numPts,numPts,columns))        
        dataSet = dataSet.mean(1)
        return dataSet
    
    #Removes all of the points in between the data points that we want to analyze
    def removeFillerPoints(self,dataSet,ppSP=1,duty=8,SPs=1):
        if duty == 0:
            return dataSet
        columns = np.shape(dataSet)[-1]
        rows = np.shape(dataSet)[0]
        period = duty+2 #One extra point for the data point and the second for the negative pulse
        rows = int(rows/period/ppSP)*period*ppSP
        dataSet = dataSet[0:rows,:]        
        dataSet = np.reshape(dataSet,(rows/period/ppSP,columns*period*ppSP))
        dataSet = dataSet[:,duty*ppSP*columns:(duty+1)*ppSP*columns]   
        dataSet = np.reshape(dataSet,(rows/(duty+2),columns))
        return dataSet
            
    def populateComboBoxes(self,dataFile):
        self.cbXAxis.clear()
        self.cbXAxis.addItems(dataFile.headers)
        self.cbYAxis.clear()        
        self.cbYAxis.addItems(dataFile.headers)
        
        self.cbAXAxis.clear()
        self.cbAXAxis.addItems(self.dataFrame.columns)
        self.cbAYAxis.clear()        
        self.cbAYAxis.addItems(self.dataFrame.columns)  
        self.cbAIndex.clear()        
        self.cbAIndex.addItems(self.dataFrame.columns)        
        
    def updateDataSet(self):
        if len(self.lwMeasSelect.selectedItems())<=0:
            self.lwMeasSelect.setCurrentRow (0)
        #Determine which measurements are selected in the list
        selectedRows = np.sort([self.lwMeasSelect.row(i) for i in self.lwMeasSelect.selectedItems()])
        #Select only the data from the desired measurements
        ppSP = self.dataFile.headerInfo['Num Points']
        setPts = self.dataFile.headerInfo['Num Setpoints']
        vdsSetPts = self.dataFile.headerInfo['Vds Setpoints']
#        df = self.dataFile.dataArray
#        duty = self.dataFile.duty
        newDataIndexes = [range(i*ppSP*setPts*vdsSetPts,(i+1)*ppSP*setPts*vdsSetPts) for i in selectedRows]
        self.dataArray_raw = self.dataFile.dataArray[np.hstack(newDataIndexes),:]      
        self.processRawData(self.dataArray_raw,ppSP,setPts)
        self.filterSweepDirection()
        self.tabChangedSlot() #Update the plot

    def updateSmoothingSlider(self,minVal,maxVal,val):
        #Must all be odd integers
        self.hsSmoothing.setMinimum((minVal-1)/2)
        self.hsSmoothing.setMaximum((maxVal-1)/2)
#        self.hsSmoothing.setSliderPosition((val-1)/2)

    def filterSweepDirection(self):
        df = self.dataFrame_orig
        #Only one sweep direction selected
        if len(self.lwSweepDirSelect.selectedItems()) <= 1:
            df = self.gradientFilter(df,self.lwSweepDirSelect.item(0).isSelected())
        df = df.sort(['Vds', 'Vgs'], ascending=[1, 1]) #Sort data to avoid discontinuities                       
        self.dataFrame = df.reindex(drop=True)
#        self.updateSmoothingSlider(minVal = 3,maxVal = 15, val = 9)        
        self.calculateMobilty()        

    #Filters dataframe according to increasing or decreasing Vgs measurements
    def gradientFilter(self,dat,forward): 
        tmp = dat.copy()
        tmp['Grad_Vgs']=np.gradient(dat['Vgs'])
        if forward:
           return dat.loc[(tmp['Grad_Vgs']>=0)]
        else:
           return dat.loc[(tmp['Grad_Vgs']<=0)]
        
    def populateMeasurementList(self,dataFile):
        self.lwMeasSelect.clear()
        for i in range(dataFile.headerInfo['Num Measurements']):
            self.lwMeasSelect.addItems(["Measurement " + str(i+1)])
        if self.numMeasurements > 1:
            self.lwMeasSelect.item(0).setSelected(True)

    def updateAnalysisInfoTables(self): #Updates table with information from analysis
        self.tblAnalysisInfo.item(self.findItemIndexes(self.tblAnalysisInfo,"Resistivity")[0],1).setText("%.4f ohm*cm"%(self.rho*1E2))
        self.tblAnalysisInfo.item(self.findItemIndexes(self.tblAnalysisInfo,"Threshold Voltage")[0],1).setText("%.2f V"%self.Vt+" +/- %.2f V"%self.VtErr)
        self.tblAnalysisInfo.item(self.findItemIndexes(self.tblAnalysisInfo,"On/Off Ratio")[0],1).setText("%.0f"%self.onOffRatio)
        self.tblAnalysisInfo.item(self.findItemIndexes(self.tblAnalysisInfo,"Peak Mobility")[0],1).setText("%.0f cm^2/Vs"%(self.peakMobility*1E4))
        self.tblAnalysisInfo.item(self.findItemIndexes(self.tblAnalysisInfo,"Electron Concentration")[0],1).setText("%.2E cm^-3"%(self.electronConcent/1E6))
        self.tblAnalysisInfo.item(self.findItemIndexes(self.tblAnalysisInfo,"Sweep Rate")[0],1).setText("%.2f mV/s"%(self.sweepRate*1E3))
        self.tblAnalysisInfo.item(self.findItemIndexes(self.tblAnalysisInfo,"Mean Free Path")[0],1).setText("%.4f nm"%(self.mfp*1E9))
        self.tblAnalysisInfo.item(self.findItemIndexes(self.tblAnalysisInfo,"Coherence Time")[0],1).setText("%.2E s"%self.cohTime)
        self.tblAnalysisInfo.item(self.findItemIndexes(self.tblAnalysisInfo,"Drift Velocity")[0],1).setText("%.0f m/s"%self.driftVelocity)         
        self.tblAnalysisInfo.item(self.findItemIndexes(self.tblAnalysisInfo,"Subthreshold Slope")[0],1).setText("%.0f mV/decade"%(self.subthresholdSlope*1E3))      
        
        self.tblAnalysisParams.item(self.findItemIndexes(self.tblAnalysisParams,"NW Length")[0],1).setText("%.0f nm"%(self.NWLength*1E9))      
        self.tblAnalysisParams.item(self.findItemIndexes(self.tblAnalysisParams,"NW Diameter")[0],1).setText("%.0f nm"%(self.NWRadius*2*1E9))     
        self.tblAnalysisParams.item(self.findItemIndexes(self.tblAnalysisParams,"Measurement Duty")[0],1).setText("%i"%(self.dataFile.duty))      
        
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
        fileName = QtGui.QFileDialog.getSaveFileName(self, 'Choose save data file name', self.dataFile.fileInfo['Filename'], filter='*.dat')
        f = open(fileName, 'w')
        f.write("# Created by MobilityCalculator - Martin Friedl\n")
        f.write("# "+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+"\n")
        if self.tabWidget.currentIndex() == 0:
            f.write('# '+self.headers[self.cbXAxis.currentIndex()]+" ("+self.units[self.cbXAxis.currentIndex()]+")"+'\t'+self.headers[self.cbYAxis.currentIndex()]+" ("+self.units[self.cbYAxis.currentIndex()]+")"+'\n')
            for x,y in zip(self.pltXData,self.pltYData):
                f.write('%E\t' % x +'%E\t' % y + '\n')            
        elif self.tabWidget.currentIndex() == 1:
            self.dataFrame.to_csv(f, sep='\t')
        f.close()            
        print "Saved to " + str(self.dataFile.filename)

    def updatePlotFromAxisSettings(self):
        #Check if the proper axes are selected to compare the Vt fit
        if (self.cbAXAxis.currentIndex() == 0) & (self.cbAYAxis.currentIndex() == 4) & (self.cbAIndex.currentIndex() == 1):
            self.ckbShowVtFit.setEnabled(True)
        else:
            self.ckbShowVtFit.setEnabled(False) 
            
        #If the Vt fit checkbox is checked, plot the Vt data first
        if (self.ckbShowVtFit.isChecked() & self.ckbShowVtFit.isEnabled()):
            self.mplwidget.myplot.hold(False)
            xaxis = self.dataFrame.columns[0]
            yaxis = self.dataFrame.columns[5]
            legend = self.dataFrame.columns[1] 
            self.updatePlot_Analysis(xaxis,yaxis,legend)  
            self.mplwidget.myplot.hold(True)            
        else:
            self.mplwidget.myplot.hold(False)
        
        xaxis = self.dataFrame.columns[self.cbAXAxis.currentIndex()]
        yaxis = self.dataFrame.columns[self.cbAYAxis.currentIndex()]
        legend = self.dataFrame.columns[self.cbAIndex.currentIndex()] 
        self.updatePlot_Analysis(xaxis,yaxis,legend)   

        self.mplwidget.myplot.hold(False)             

    def updatePlot_Analysis(self,xaxis,yaxis,legend):
        try:
            plotFrame = pd.pivot_table(self.dataFrame,values=yaxis, index=xaxis, columns=legend)
        except:
            return
        #if data set has more legend entries than is allowed
        if plotFrame.columns.size > self.sbLegendEntries.value():
            #Create equally spaced indexes that will be put in the legend
            legendIndexes = [int(np.round(i)) for i in np.linspace(0,plotFrame.columns.size-1,self.sbLegendEntries.value())]
            #Select only this data to then plot            
            plotFrame = plotFrame.iloc[:,legendIndexes]

        self.pltCfg_Analysis['xlabel'] = xaxis
        self.pltCfg_Analysis['ylabel'] = yaxis      
        self.pltCfg_Analysis['title'] = self.dataFile.fileInfo['Wafer']+" "+self.dataFile.fileInfo['Device ID']+" Exp:"+self.dataFile.fileInfo['Measurement']
        self.plotDataFrame(plotFrame,self.pltCfg_Analysis)                          
        
    def plotDataFrame(self,dataFrame,pltCfg):
        
        holdPlot = self.mplwidget.myplot.ishold()
        if not(holdPlot): #If hold plot is off, clear the plot before plotting new one
            self.mplwidget.myplot.cla()  
        self.mplwidget.myplot.hold(True) #Turn on hold to plot multiple data sets
        dataFrame.plot(ax=self.mplwidget.myplot)
        self.mplwidget.myplot.plot()
        
        self.mplwidget.myplot.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
        self.mplwidget.myplot.ticklabel_format(style='sci', axis='x', scilimits=(-2,2))  
        self.mplwidget.myplot.set_xlabel(pltCfg['xlabel'])
        self.mplwidget.myplot.set_ylabel(pltCfg['ylabel'])
        self.mplwidget.myplot.set_title(pltCfg['title'])        
        self.lineOrMarkerStyleChangedSlot()#Update line settings to match user selection      
        self.mplwidget.figure.canvas.draw()
        
        #Reset the toolbar home button for the new plot
        self.mplwidget.mpl_toolbar._views.clear()
        self.mplwidget.mpl_toolbar._positions.clear()        
        
        if not(holdPlot): #If hold plot is off, return plot back to original (not hold) setting
            self.mplwidget.myplot.hold(False)        
            
    def plotData(self,xpoints,ypoints,pltCfg_RawData):
        self.mplwidget.myplot.cla()
        self.mplwidget.myplot.plot(xpoints,ypoints,marker=pltCfg_RawData['marker'],linestyle=pltCfg_RawData['linestyle'], color=pltCfg_RawData['linecolor'])
        self.mplwidget.myplot.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
        self.mplwidget.myplot.ticklabel_format(style='sci', axis='x', scilimits=(-2,2))  
        self.mplwidget.myplot.set_xlabel(pltCfg_RawData['xlabel'])
        self.mplwidget.myplot.set_ylabel(pltCfg_RawData['ylabel'])
        self.mplwidget.myplot.set_title(pltCfg_RawData['title'])
        self.mplwidget.myplot.grid()
        self.mplwidget.figure.canvas.draw()          
        
    def updatePlot_RawData(self):
        try:
            if self.ckbDriftCorr.isChecked():#If doing drift correction, automatically filter and average
                self.pltXData = self.dataArray_avgFilt_driftCorr[:,self.cbXAxis.currentIndex()]
                self.pltYData = self.dataArray_avgFilt_driftCorr[:,self.cbYAxis.currentIndex()]
            else:
                if self.ckbAvgData.isChecked() and not(self.ckbRmvFiller.isChecked()):
                    self.pltXData = self.dataArray_avg[:,self.cbXAxis.currentIndex()]
                    self.pltYData = self.dataArray_avg[:,self.cbYAxis.currentIndex()]            
                elif not(self.ckbAvgData.isChecked()) and self.ckbRmvFiller.isChecked():
                    self.pltXData = self.dataArray_filt[:,self.cbXAxis.currentIndex()]
                    self.pltYData = self.dataArray_filt[:,self.cbYAxis.currentIndex()] 
                elif self.ckbAvgData.isChecked() and self.ckbRmvFiller.isChecked():
                    self.pltXData = self.dataArray_avgFilt[:,self.cbXAxis.currentIndex()]
                    self.pltYData = self.dataArray_avgFilt[:,self.cbYAxis.currentIndex()] 
                else:
                    self.pltXData = self.dataArray_raw[:,self.cbXAxis.currentIndex()]
                    self.pltYData = self.dataArray_raw[:,self.cbYAxis.currentIndex()]       
        except:
            return
            
        self.mplwidget.myplot.cla()
        marker = self.markerDict[self.cbMarkerStyle.currentIndex()]
        linestyle = self.lineStyleDict[self.cbLineStyle.currentIndex()]
        colour = self.lineColourDict[self.cbPlotColour.currentIndex()]
        self.mplwidget.myplot.plot(self.pltXData,self.pltYData,marker=marker, linestyle=linestyle, color=colour)
        self.mplwidget.myplot.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
        self.mplwidget.myplot.ticklabel_format(style='sci', axis='x', scilimits=(-2,2))
        self.mplwidget.myplot.set_xlabel(self.dataFile.headers[self.cbXAxis.currentIndex()]+" ("+self.dataFile.units[self.cbXAxis.currentIndex()]+")")
        self.mplwidget.myplot.set_ylabel(self.dataFile.headers[self.cbYAxis.currentIndex()]+" ("+self.dataFile.units[self.cbYAxis.currentIndex()]+")")
        self.mplwidget.myplot.set_title(self.dataFile.fileInfo['Wafer']+" "+self.dataFile.fileInfo['Device ID']+" Exp:"+self.dataFile.fileInfo['Measurement'])
        self.mplwidget.myplot.grid()
        self.mplwidget.figure.canvas.draw()

        #Reset the toolbar home button for the new plot
        self.mplwidget.mpl_toolbar._views.clear()
        self.mplwidget.mpl_toolbar._positions.clear()        
        
if __name__ == '__main__':
    # Create the GUI application
    qApp = QtGui.QApplication(sys.argv)
    window = MyWindow()
    window.show() 
#    window.loadData("D:/Dropbox/LMSC/Programming/ACMobility/MF-04-IA_04-23B_96-26-0_dark_MobilityPulsedAC-AD.lvm")
#    window.loadData("D:/Dropbox/LMSC/Programming/ACMobility/MF-04-IA_04-23B_92-26-1_dark_MobilityPulsedAC-AD.lvm")
#    window.loadData("./MF-04-IA_04-23B_88-16-2_dark_MobilityPulsedAC-AD.lvm")

#    window.loadData("./MF-04-IA_04-23B_74-96-2_dark_MobilityPulsedAC-AD.lvm")

#    window.loadData("./MF-08-ZA_00-00_84-91-2_dark_MobilityPulsedAC-BC.lvm")

#    window.loadData("./MF-06-IA_04-23B_20-97-1_dark_MobilityPulsedAC-AD.lvm")
    window.loadData("V:/2015-04-09 - MF-06-IA - Mobility Measurement (Python)/MF-06-IA_00-00_41-27-6_dark_LinearMobilityAD.dat")

#    testFile = dataFile("D:/Dropbox/LMSC/Programming/ACMobility/MF-04-IA_04-23B_82-26-1_dark_MobilityPulsedAC-AD-VdsFirst.lvm")    


    window.tabWidget.setCurrentIndex(1)
    window.cbAXAxis.setCurrentIndex(0)
    window.cbAYAxis.setCurrentIndex(4)
    window.cbAIndex.setCurrentIndex(1)

    # start the Qt main loop execution, exiting from this script
    # with the same return code of Qt application
    sys.exit(qApp.exec_()) 
    
    