# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 09:18:21 2015

@author: Martin Friedl
"""
import re
# Python Qt4 bindings for GUI objects
from PyQt4 import QtCore,QtGui, uic
# Numpy functions for image creation
import numpy as np
# Pandas for manipulating datasets
import pandas as pd

class dataFile(object):
    def __init__(self, filename = None):
        #Filename Information
        tmp = { 'Filename' : filename,\
                'Wafer' : None,\
                'NW Type' : None,\
                'Device ID' : None,\
                'Measurement' : None,\
                'Lighting' : None,\
                'Test Type' : None}
        self.fileInfo = pd.Series(tmp)
        
        #Header File Information
        tmp = { 'Delay Init' : None,\
                'Delay After Set' : None,\
                'Delay After Read' : None,\
                'Num Points' : None,\
                'Num Measurements' : None,\
                'Num Setpoints' : None,\
                'Vds Setpoints' : None}
        self.headerInfo = pd.Series(tmp)        
        
        #DataSet Initialization
        self.headers = None
        self.units = None        
        self.duty = 8       
        self.dataArray = None     
        
        if not(filename is None):
            self.load(filename)
        
    def load(self,fileName = None):
        if fileName is None:
            fileName = QtGui.QFileDialog.getOpenFileName(self, 'Choose data file to open', self.fileInfo['Filename'], filter='Measurement Files (*.dat *.lvm);;All Files (*.*)')
        if fileName:
            self.filename = fileName
            self.parseFilename()
            if "linear" in self.fileInfo['Test Type'].lower(): #If linear sweep, duty is 0
                self.duty = 0
            f = open(self.filename, 'r')
            self.readFileHeader(f)
            f.seek(0) #Reset cursor to beginning of the file
            self.dataArray = self.importData(f) #Original data array that never changes
            f.close()          
    
    def parseFilename(self): #Extract information from the data filename
        parsed = str(self.filename.split("/")[-1]).split(".")[0].split("_")
        self.fileInfo['Wafer'] = parsed[0]
        self.fileInfo['NW Type'] = parsed[1]
        self.fileInfo['Device ID'] = parsed[2].split("-")[0]+"-"+parsed[2].split("-")[1]
        self.fileInfo['Measurement'] = parsed[2].split("-")[-1]
        self.fileInfo['Lighting'] = parsed[3]
        self.fileInfo['Test Type'] = parsed[4]
        
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
                        self.headerInfo['Delay Init'] = value
                    elif "delay(after set)" in text:
                        self.headerInfo['Delay After Set'] = value
                    elif "delay(bef.read)" in text:
                        self.headerInfo['Delay After Read'] = value                    
                    elif "np. of points" in text:
                        self.headerInfo['Num Points'] = value                   
                    elif "np. of measurements" in text:
                        self.headerInfo['Num Measurements'] = value                   
                    elif "np. of setpoints" in text:
                        self.headerInfo['Num Setpoints'] = value             
                    elif float(text) == value:
                        self.headerInfo['Vds Setpoints'] = value
                    else:
                        pass    
        
    def importData(self,f):
        self.headers = None
        self.units = None
        dataArray = None
        lineDataList = []
        for line in f: 
            if (line[0] == '#')or(line == '\n'): continue
            splitLine = line[:-1].split('\t')
            splitLine = filter(bool, splitLine) #Remove empty entries in the list
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
        #Determine number of complete measurements   
        datapoints = np.shape(dataArray)[0]        
        self.headerInfo['Num Measurements'] = int(datapoints/1.0/self.headerInfo['Num Points']/self.headerInfo['Num Setpoints']/self.headerInfo['Vds Setpoints'])
        #Truncate data to exclude incomplete measurements
        dataArray = dataArray[0:self.headerInfo['Num Measurements']*self.headerInfo['Num Points']*self.headerInfo['Num Setpoints']*self.headerInfo['Vds Setpoints'],:]        
        return dataArray     
        
    def updateDeviceInfoTable(self,tbl): #Updates table with information from filename of data file
        self.updateTable(tbl,self.fileInfo)
        
    def updateHeaderInfoTable(self,tbl): #Updates table with information from header of data file
        self.updateTable(tbl,self.headerInfo)
    
    def updateTable(self,tbl,infoSeries):
        if not(type(tbl) is QtGui.QTableWidget):
            raise Exception('Expected type QTableWidget!')
        tbl.clear()
        tbl.setRowCount(0)
        tbl.setHorizontalHeaderLabels(['Parameter','Value'])
        
        for i in range(len(infoSeries)):
            tbl.insertRow(i)            
            item1 = QtGui.QTableWidgetItem(str(infoSeries.index[i]))
            tbl.setItem(i, 0, item1)            
            item2 = QtGui.QTableWidgetItem(str(infoSeries.loc[infoSeries.index[i]]))
            tbl.setItem(i, 1, item2)        
        
if __name__ == '__main__':
    testFile = dataFile("D:/Dropbox/LMSC/Programming/ACMobility/MF-04-IA_04-23B_88-16-2_dark_MobilityPulsedAC-AD.lvm")    
    
    print testFile.dataArray
    print testFile.fileInfo
    bsdf = testFile.fileInfo
    sdfa = testFile.headerInfo