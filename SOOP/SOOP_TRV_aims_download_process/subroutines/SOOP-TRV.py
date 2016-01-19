#/usr/bin/env python
# Download SOOP TRV data from AIMS Web Service
# The script reads an XML file provided by AIMS. The script then looks at which new
# channel is available to download, and compare this list with a pickle file (a python
# way to store python variables) containing what has already been downloaded.
# Some modifications on the files have to be done in order to be CF and IMOS compliant
# The files are stored in data_wip_path as defined by confix.txt
#
# author Laurent Besnard, laurent.besnard@utas.edu.au

import urllib2
import urllib
import xml.etree.ElementTree as ET
import tempfile
import time
import zipfile
import logging
import pickle
import os
import subprocess, shlex
import shutil
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from destPath import *
from time import gmtime, strftime
from datetime import date

#####################
# Logging Functions #
#####################
# start logging using logging python library
# output:
#     logger - similar to a file handler
def loggingAims():
    wipPath                   = os.environ.get('data_wip_path')

    logging.basicConfig(level = logging.INFO)
    logger                    = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # create a file handler
    handler                   = logging.FileHandler(os.path.join(wipPath, 'soop_trv.log'))
    handler.setLevel(logging.INFO)

    # create a logging format
    formatter                 = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(handler)
    return logger

# close logging
# input:
#     logger : logging handler generated by loggingAims()
def closeLogger(logger):
    #closes the handlers of the specified logger only
    x = list(logger.handlers)
    for i in x:
        logger.removeHandler(i)
        i.flush()
        i.close()

####################
# Pickle Functions #
####################
# returns the pickle filepath according to the QC level being processed
# input:
#     levelQc(int) : 0 or 1
# output:
#     picleQcFile(str) : pickle file path
def _pickleFileName(levelQc):
    wipPath = os.environ.get('data_wip_path')

    if levelQc == 0:
        pickleQcFile = os.path.join(wipPath, 'soop_trv_qc0.pickle')
    elif levelQc == 1:
        pickleQcFile = os.path.join(wipPath, 'soop_trv_qc1.pickle')

    return pickleQcFile

# if channelId has been successfuly processed, we write about it in a pickle file
# the the channelId doesn't get reprocessed
# input:
#    channelId(str)      : channelId to save information
#    aimsXmlInfo(tupple) : generated by parserAimsXml
#    levelQc(int)        : 0 or 1
def saveChannelInfo(channelId, aimsXmlInfo, levelQc):
    pickleFile = _pickleFileName(levelQc)

    # condition in case the pickle file already exists or not. In the first case,
    # aimsXmlInfo comes from the pickle, file, otherwise comes from the function arg
    if os.path.isfile(pickleFile):
        with open(pickleFile, 'rb') as pRead:
            aimsXmlInfo        = pickle.load(pRead)

        channelIdIndex     = aimsXmlInfo[0].index(channelId) # value important to write lastDownloadedDate to correct index
        lastDownloadedDate = aimsXmlInfo[-1]
    else:
        lastDownloadedDate = [None] * len(aimsXmlInfo[0]) # initialise array
        channelIdIndex     = aimsXmlInfo[0].index(channelId)

    channelIdInfo                      = getChannelInfo(channelId, aimsXmlInfo)
    lastDownloadedDate[channelIdIndex] = channelIdInfo[2] # fromDate
    pickleDb                           = aimsXmlInfo[0:-1] + (lastDownloadedDate,) # add to tupple

    with open(pickleFile, 'wb') as pWrite:
        pickle.dump(pickleDb, pWrite)

def hasChannelAlreadyBeenDownloaded(channelId, levelQc):
    pickleFile = _pickleFileName(levelQc) # different pickle per QC
    if os.path.isfile(pickleFile):
        with open(pickleFile, 'rb') as pRead:
            pickleDb = pickle.load(pRead)

        if channelId in pickleDb[0]: # check the channel is in the pickle file
            channelIdIndex = pickleDb[0].index(channelId)

            if pickleDb[-1][channelIdIndex] is not None: # check the last downloadedDate field
                return True
            else:
                return False
        else:
            return False

    else:
        return False

#############
# UNIT TEST #
#############
import hashlib
def md5(fname):
    hash = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash.update(chunk)
    return hash.hexdigest()

# Check that a the AIMS system or this script hasn't been modified. This function checks that a downloaded file still has the same ch5.
def unitTestSoopTrv():
    channelId               = '48594'
    fromDate                = '2014-08-30T14:51:34Z'
    thruDate                = '2014-09-01T14:50:51Z'
    levelQc                 = '1'
    xmlUrl                  = 'http://data.aims.gov.au/gbroosdata/services/rss/netcdf/level' + str(levelQc) +'/100'

    md5Value                = 'da0b8c836fe8ba7d167fcac05efd8ea8'
    aimsXmlInfo             = parseAimsSoopTrvXml(xmlUrl)
    channelIdInfo           = getChannelInfo(channelId, aimsXmlInfo)
    netcdfTmpFilePath       = downloadChannel(channelId, fromDate, thruDate, levelQc)
    modifySoopTrvNetcdf(netcdfTmpFilePath, channelIdInfo)

    # force values of attributes which change all the time
    netcdfFileObj              = Dataset(netcdfTmpFilePath, 'a', format='NETCDF4')
    netcdfFileObj.date_created = "2015-11-25T05:03:25Z"
    netcdfFileObj.history      = 'unit test only'
    netcdfFileObj.close()
    booleanReturn              = md5(netcdfTmpFilePath) == md5Value
    os.remove(netcdfTmpFilePath)

    return booleanReturn

######################
# XML Info Functions #
######################
def parseAimsSoopTrvXml(xmlUrl):
    logger.info('parse SOOP-TRV xml : ' + xmlUrl)
    response      = urllib2.urlopen(xmlUrl)
    html          = response.read()
    root          = ET.fromstring(html)

    nextItemExist = True
    nItem         = 3 # start number for AIMS xml file

    title         = []
    link          = []
    metadataUuid  = []
    uom           = []
    fromDate      = []
    thruDate      = []
    platformName  = []
    siteName      = []
    channelId     = []
    parameter     = []
    parameterType = []
    tripId        = []

    while nextItemExist:
        title        .append(root[0][nItem][0].text)
        link         .append(root[0][nItem][1].text)
        metadataUuid .append(root[0][nItem][6].text)
        uom          .append(root[0][nItem][7].text)
        fromDate     .append(root[0][nItem][8].text)
        thruDate     .append(root[0][nItem][9].text)
        platformName .append(root[0][nItem][10].text)
        siteName     .append(root[0][nItem][11].text)
        channelId    .append(root[0][nItem][12].text)
        parameter    .append(root[0][nItem][13].text)
        parameterType.append(root[0][nItem][14].text)

        # in case there is no trip id defined by AIMS, we create a fake one
        try:
            tripId       .append(root[0][nItem][15].text)
        except IndexError:
            dateObject = time.strptime(root[0][nItem][8].text,"%Y-%m-%dT%H:%M:%SZ")
            tripIdFake =  str(dateObject.tm_year) + str(dateObject.tm_mon).zfill(2) + str(dateObject.tm_mday).zfill(2)
            tripId.append(tripIdFake)

        nItem += 1
        # test if next item in XML file exists
        try:
            root[0][nItem +1]
            nextItemExist = True
        except IndexError:
            nextItemExist = False

    return channelId, fromDate, thruDate, metadataUuid, uom, platformName, siteName, parameter, parameterType, tripId

# returns all the informations found in the parsed AIMS xml file for one channel only
# input:
#    channelId(int) : channelId to return info from
#    aimsXmlInfo    : generated by parseAimsSoopTrvXml
# output:
#    similar as parseAimsSoopTrvXml
def getChannelInfo(channelId, aimsXmlInfo):
    channelIdIndex = aimsXmlInfo[0].index(channelId)
    fromDate       = aimsXmlInfo[1][channelIdIndex]
    thruDate       = aimsXmlInfo[2][channelIdIndex]
    metadataUuid   = aimsXmlInfo[3][channelIdIndex]
    uom            = aimsXmlInfo[4][channelIdIndex]
    platformName   = aimsXmlInfo[5][channelIdIndex]
    siteName       = aimsXmlInfo[6][channelIdIndex]
    parameter      = aimsXmlInfo[7][channelIdIndex]
    parameterType  = aimsXmlInfo[8][channelIdIndex]
    tripId         = aimsXmlInfo[9][channelIdIndex]

    return channelId, fromDate, thruDate, metadataUuid, uom, platformName, siteName, parameter, parameterType, tripId

##########################################
# Channel Process/Download/Mod Functions #
##########################################
# generated the data link to download, and extract the zip file into a temp file
# input:
#    channelId(str) : channelId to download
#    fromDate(str)  : str containing the first time to start the download from written in this format 2009-04-21T10:43:54Z
#    thruDate(str)  : same as above but for the last date
#    levelQc(int)   : 0 or 1
def downloadChannel(channelId, fromDate, thruDate, levelQc):
    urlDataDownload = 'http://data.aims.gov.au/gbroosdata/services/data/rtds/' + channelId + '/' + 'level' + str(levelQc) + '/raw/raw/' + fromDate + '/' + thruDate + '/netcdf/2'
    tmpZipFile      = tempfile.mkstemp()
    tmpZipFileName  = tmpZipFile[1]
    urllib.urlretrieve(urlDataDownload, tmpZipFileName)

    netcdfTmpPath   = tempfile.mkdtemp()
    zip             = zipfile.ZipFile(tmpZipFileName)

    for name in zip.namelist():
        zip.extract(name, netcdfTmpPath)
        netcdfFilePath = os.path.join(netcdfTmpPath, name)

    zip.close()
    os.remove(tmpZipFileName)
    logger.info("    %s downloaded successfuly" % urlDataDownload)
    return netcdfFilePath

# Check if the unzipped file is a 'NO_DATA_FOUND' file instead of a NetCDF file
# this behaviour is correct for FAIMMS and NRS, as it means no data for the selected
# time period. However it doesn't make sense for SOOP TRV, This should be treated as
# an error
def _isNoDataFound(netcdfFilePath):
    return os.path.basename(netcdfFilePath) == 'NO_DATA_FOUND'

# Rename global attribute from netcdf4 dataset object
#  object             = Dataset(netcdfFile, 'a', format='NETCDF4')
#  old_attribute_name = current gatt name to modify
#  new_attribute_name = new gatt name
def renameNetcdfAttribute(object_, old_attribute_name, new_attribute_name):
    setattr(object_, new_attribute_name, getattr(object_, old_attribute_name))
    delattr(object_, old_attribute_name)

# some channels have only _Fillvalues in their main variable. This is not correct and need
# to be tested
def hasMainVarOnlyFillValue(netcdfFilePath):
    mainVar       = getMainSoopTrvVar(netcdfFilePath)
    netcdfFileObj = Dataset(netcdfFilePath, 'a', format='NETCDF4')
    varObj        = netcdfFileObj.variables[mainVar]
    varValues     = varObj[:]

    # if no fill value in variable, no mask attribute
    if hasattr(varValues,'mask'):
        return all(varValues.mask)
    else:
        return False

# convert a CF time into an IMOS one forced to be 'days since 1950-01-01 00:00:00'
# the variable HAS to be 'TIME'
def convertTimeCftoImos(netcdfFilePath):
    try:
        netcdfFileObj = Dataset(netcdfFilePath, 'a', format='NETCDF4')
        time          = netcdfFileObj.variables['TIME']
        dtime         = num2date(time[:], time.units, time.calendar)   # this gives an array of datetime objects
        time.units    = 'days since 1950-01-01 00:00:00 UTC'
        time[:]       = date2num(dtime, time.units, time.calendar) # conversion to IMOS recommended time
        netcdfFileObj.close()
        return True
    except:
        return False

# Modify the downloaded NetCDF file so it passes both CF and IMOS checker
# input:
#    netcdfFilePath(str)    : path of netcdf file to modify
#    channelIdIndex(tupple) : information from xml for the channel
def modifySoopTrvNetcdf(netcdfFilePath, channelIdInfo):
    netcdfFileObj = Dataset(netcdfFilePath, 'a', format='NETCDF4')
    shipCode      = netcdfFileObj.platform_code

    if shipCode   == 'VNCF':
        vesselName   = 'Cape-Ferguson'
    elif shipCode == 'VMQ9273':
        vesselName   = 'Solander'
    else:
        logger.error('   UNKNOWN SHIP - channel ' + str(channelIdInfo[0]))
        netcdfFileObj.close()
        return False

    # add gatts to NetCDF
    netcdfFileObj.aims_channel_id = int(channelIdInfo[0])
    netcdfFileObj.cdm_data_type   = 'Trajectory'
    netcdfFileObj.vessel_name     = vesselName
    netcdfFileObj.trip_id         = int(channelIdInfo[9])

    if not (channelIdInfo[3] == 'Not Available'):
        netcdfFileObj.metadata_uuid = channelIdInfo[3]

    # add CF gatts, values stored in lib/netcdf/netcdf-cf-imos-compliance.sh
    netcdfFileObj.Conventions            = os.environ.get('CONVENTIONS')
    netcdfFileObj.data_centre_email      = os.environ.get('DATA_CENTRE_EMAIL')
    netcdfFileObj.data_centre            = os.environ.get('DATA_CENTRE')
    netcdfFileObj.project                = os.environ.get('PROJECT')
    netcdfFileObj.acknowledgement        = os.environ.get('ACKNOWLEDGEMENT')
    netcdfFileObj.distribution_statement = os.environ.get('DISTRIBUTION_STATEMENT')

    netcdfFileObj.cdm_data_type       = "Trajectory"
    netcdfFileObj.date_created        = strftime("%Y-%m-%dT%H:%M:%SZ", gmtime())
    netcdfFileObj.quality_control_set = 1
    imosQcConvention                  = 'IMOS standard set using the IODE flags'
    netcdfFileObj.author              = 'laurent besnard'
    netcdfFileObj.author_email        = 'laurent.besnard@utas.edu.au'

    renameNetcdfAttribute(netcdfFileObj, 'geospatial_LAT_max', 'geospatial_lat_max')
    renameNetcdfAttribute(netcdfFileObj, 'geospatial_LAT_min', 'geospatial_lat_min')
    renameNetcdfAttribute(netcdfFileObj, 'geospatial_LON_max', 'geospatial_lon_max')
    renameNetcdfAttribute(netcdfFileObj, 'geospatial_LON_min', 'geospatial_lon_min')

    # variables modifications
    time           = netcdfFileObj.variables['time']
    time.calendar  = 'gregorian'
    time.axis      = 'T'
    time.valid_min = 0.0
    time.valid_max = 9999999999.0
    netcdfFileObj.renameVariable('time','TIME')
    netcdfFileObj.renameDimension('time','TIME')

    # depth
    depth                 = netcdfFileObj.variables['depth']
    depth.positive        = 'down'
    depth.axis            = 'Z'
    depth.reference_datum = 'graphical coordinates, WGS84 projection'
    depth.valid_max       = 30.0
    depth.valid_min       = -10.0
    netcdfFileObj.renameVariable('depth','DEPTH')

    # latitude longitude
    latitude                     = netcdfFileObj.variables['LATITUDE']
    latitude.axis                = 'Y'
    latitude.valid_min           = -90.0
    latitude.valid_max           = 90.0
    latitude.reference_datum     = 'geographical coordinates, WGS84 projection'
    latitude.standard_name       = 'latitude'
    latitude.long_name           = 'latitude'
    latitude.ancillary_variables = 'LATITUDE_quality_control'

    longitude                     = netcdfFileObj.variables['LONGITUDE']
    longitude.axis                = 'X'
    longitude.valid_min           = -180.0
    longitude.valid_max           = 180.0
    longitude.reference_datum     = 'geographical coordinates, WGS84 projection'
    longitude.standard_name       = 'longitude'
    longitude.long_name           = 'longitude'
    longitude.ancillary_variables = 'LONGITUDE_quality_control'

    latitude_qc                              = netcdfFileObj.variables['LATITUDE_quality_control']
    latitude_qc.standard_name                = 'latitude status_flag'
    latitude_qc.quality_control_conventions  = imosQcConvention
    longitude_qc                             = netcdfFileObj.variables['LONGITUDE_quality_control']
    longitude_qc.standard_name               = 'longitude status_flag'
    longitude_qc.quality_control_conventions = imosQcConvention

    if 'fluorescence' in netcdfFileObj.variables.keys():
        netcdfFileObj.renameVariable('fluorescence','CPHL')

        cphl             = netcdfFileObj.variables['CPHL']
        cphl.long_name   = 'mass_concentration_of_inferred_chlorophyll_from_relative_fluorescence_units_in_sea_waterconcentration_of_chlorophyll_in_sea_water'
        cphl.coordinates = "TIME LATITUDE LONGITUDE DEPTH"
        delattr(cphl, 'standard_name')

        netcdfFileObj.renameVariable('fluorescence_quality_control','CPHL_quality_control')
        netcdfFileObj.variables['CPHL_quality_control'].long_name = 'mass_concentration_of_inferred_chlorophyll_from_relative_fluorescence_units_in_sea_waterconcentration_of_chlorophyll_in_sea_water status_flag'
        delattr(netcdfFileObj.variables['CPHL_quality_control'], 'standard_name')

        cphl.ancillary_variables                                                    = 'CPHL_quality_control'
        netcdfFileObj.variables['CPHL_quality_control'].quality_control_conventions = imosQcConvention

    elif 'Seawater_Intake_Temperature' in netcdfFileObj.variables.keys():
        temp                                                                     = netcdfFileObj.variables['Seawater_Intake_Temperature']
        temp.units                                                               = 'Celsius'
        temp.coordinates                                                         = "TIME LATITUDE LONGITUDE DEPTH"
        netcdfFileObj.renameVariable('Seawater_Intake_Temperature','TEMP')
        netcdfFileObj.renameVariable('Seawater_Intake_Temperature_quality_control','TEMP_quality_control')
        temp.ancillary_variables                                                 = 'TEMP_quality_control'
        netcdfFileObj.variables['TEMP_quality_control'].quality_control_conventions = imosQcConvention

    elif 'PSAL' in netcdfFileObj.variables.keys():
        netcdfFileObj.variables['PSAL'].units                                       = '1e-3'
        netcdfFileObj.variables['PSAL'].coordinates                                 = "TIME LATITUDE LONGITUDE DEPTH"
        netcdfFileObj.variables['PSAL_quality_control'].quality_control_conventions = imosQcConvention

    elif 'TURB' in netcdfFileObj.variables.keys():
        turb                                                                        = netcdfFileObj.variables['TURB']
        turb.units                                                                  = '1'
        turb.standard_name                                                          = 'sea_water_turbidity'
        turb.coordinates                                                            = "TIME LATITUDE LONGITUDE DEPTH"
        netcdfFileObj.variables['TURB_quality_control'].standard_name               = 'sea_water_turbidity status_flag'
        netcdfFileObj.variables['TURB_quality_control'].quality_control_conventions = imosQcConvention

    netcdfFileObj.close()

    if not convertTimeCftoImos(netcdfFilePath):
        return False

    _removeDimensionFromNetcdf(netcdfFilePath) # last modification to do !

    logger.info('    Modified NetCDF file to pass CF and IMOS Compliance Checker')
    return True

# DIRTY, calling bash. need to write in Python, or part of the NetCDF4 module
def _removeDimensionFromNetcdf(netcdfFilePath):
    # need to remove the 'single' dimension name from DEPTH. Unfortunately can't seem to find a way to do it easily with netCDF4 module
    subprocess.check_call(['ncwa','-O','-a','single',netcdfFilePath,netcdfFilePath])

# Some files had in the past bad latitude/longitude tracks. This is a
# really easy way to check this
# netcdfFilePath9str) : path of the netcdf file to check
def _isLatLonValuesOutsideBoundaries(netcdfFilePath):
    netcdfFileObj = Dataset(netcdfFilePath, 'a', format='NETCDF4')
    lat           = netcdfFileObj.variables['LATITUDE'][:]
    lon           = netcdfFileObj.variables['LONGITUDE'][:]
    netcdfFileObj.close()

    return any(lat > 0) or any(lat < -50) or any(lon > 180) or any(lon < 0)

def moveNetcdfToIncomingDir(netcdfPath):
    incomingDir     = os.environ.get('INCOMING_DIR')
    soopIncomingDir = os.path.join(incomingDir,'SOOP/TRV',os.path.basename(netcdfPath))

    shutil.copy(netcdfPath, soopIncomingDir) # WARNING, shutil.move creates a wrong incron event
    os.remove(netcdfPath)
    logger.info('    NetCDF moved to INCOMING DIRECTORY')

# Downloads all the data available for one channelId and moves the file to a wipPath dir
# channelId(str)
# aimsXmlInfo(tuple)
# levelQc(int)
def processChannel(channelId, aimsXmlInfo, levelQc):
    channelIdInfo = getChannelInfo(channelId,aimsXmlInfo)
    if not hasChannelAlreadyBeenDownloaded(channelId,levelQc):
        logger.info('>> QC' + str(levelQc) + ' - Processing channel ' + str(channelId))
        fromDate           = channelIdInfo[1]
        thruDate           = channelIdInfo[2]
        netcdfTmpFilePath  = downloadChannel(channelId,fromDate,thruDate,levelQc)

        if _isNoDataFound(netcdfTmpFilePath):
            logger.error('   Channel ' + str(channelId) + ' - NO_DATA_FOUND file in Zip file - CONTACT AIMS')
            os.remove(netcdfTmpFilePath)
            os.rmdir(os.path.dirname(netcdfTmpFilePath))
            return False

        if not modifySoopTrvNetcdf(netcdfTmpFilePath,channelIdInfo):
            os.remove(netcdfTmpFilePath)
            os.rmdir(os.path.dirname(netcdfTmpFilePath))
            return False

        if hasMainVarOnlyFillValue(netcdfTmpFilePath):
            logger.error('   Channel ' + str(channelId) + ' - _FillValues only in main variable - CONTACT AIMS')
            os.remove(netcdfTmpFilePath)
            os.rmdir(os.path.dirname(netcdfTmpFilePath))
            return False

        if _isLatLonValuesOutsideBoundaries(netcdfTmpFilePath):
            logger.error('   Channel ' + str(channelId) + ' - Lat/Lon values outside of boundaries - CONTACT AIMS')
            os.remove(netcdfTmpFilePath)
            os.rmdir(os.path.dirname(netcdfTmpFilePath))
            return False

        moveNetcdfToIncomingDir(netcdfTmpFilePath)
        return True

    else:
        logger.info('QC' + str(levelQc) + ' - Channel ' + str(channelId) + ' already processed')
        return False

# Downloads all channels for a QC level
# levelQc(int) : 0 or 1
def processQcLevel(levelQc):
    logger.info('Process SOOP-TRV download from AIMS web service - QC level ' + str(levelQc))
    xmlUrl      = 'http://data.aims.gov.au/gbroosdata/services/rss/netcdf/level' + str(levelQc) +'/100'
    aimsXmlInfo = parseAimsSoopTrvXml(xmlUrl)

    for channelId in aimsXmlInfo[0]:
        try:
            isChannelProcessed = processChannel(channelId,aimsXmlInfo,levelQc)
            if isChannelProcessed:
                saveChannelInfo(channelId,aimsXmlInfo,levelQc)
        except:
            logger.error('   Channel ' + str(channelId) +  ' QC' + str(levelQc) + ' - Failed, unknown reason - manual debug required')


if __name__=='__main__':

    wipPath = os.environ.get('data_wip_path')
    if not wipPath:
        logger.error('data_wip_path from config.txt is empty')
        exit

    if not os.path.exists(wipPath):
        os.makedirs(wipPath)

    logger  = loggingAims()
    unitTestResult= unitTestSoopTrv()

    if unitTestResult == True:
        logger.info('Unit Test passed')
        processQcLevel(1) # no need to process level 0 for SOOP TRV
    else:
        logger.error('Unit Test failed - file changed by AIMS or emII')

    closeLogger(logger)
    exit(0)
