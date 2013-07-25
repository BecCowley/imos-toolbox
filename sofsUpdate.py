#! /usr/bin/env python
#
# Get latest SOFS files from BOM ftp site and sort them into the
# appropriate directories on opendap.


from IMOSfile.IMOSfilename import parseFilename
import re
from datetime import datetime
import os
import argparse


def dataCategory(info):
    """
    Return the appropriate data category for a file with the given
    attributes (info as returned by parseFilename()).
    """

    if info['data_code'] == 'CMST':
        return 'Surface_properties'

    if info['data_code'] == 'FMT':
        return 'Surface_fluxes'

    if info['data_code'] == 'W':
        return 'Surface_waves'

    if info['data_code'] == 'V':
        return 'Sub-surface_currents'

    if info['data_code'] in ['PT', 'CPT']:
        return 'Sub-surface_temperature_pressure_conductivity'

    return '???'
    

def destPath(info, basePath=''):
    """
    Return the pubplic directory path for a file with the given info
    (as returned by parseFilename()).
    """
    if not info.has_key('data_category'):
        info['data_category'] = dataCategory(info)

    if info['data_category'] == '???':
        return ''

    path = os.path.join(basePath, info['data_category'])
    
    if type(info['end_time']) is datetime:
        # whole-deployment delayed-mode file, no further sub-directories
        return path

    if info['product_code'] != '1-min-avg':
        # i.e. not a daily delayed-mode product
        path = os.path.join(path, 'Real-time')

    path = os.path.join(path, info['start_time'].strftime('%Y') + '_daily' )

    return path


def updateFile(source, destDir, log=None, dry_run=False):
    """
    Synchronise source file to destDir, copying the file only if it
    doesn't exist at destDir or if source has been modified more
    rencently. If destDir is not an existing directory, create it.  If
    log is set to an open file object, log the sync command to it.  If
    dry_run is True, don't more anything, just log actions.

    Return a single character indicating the result: 
    'A' - added new file; 
    'U' - updated existing file; 
    'E' - existing file, no update;
    'F' - update failed.
    """
    # make sure source is a file
    if not os.path.isfile(source):
        print 'WARNING (updateFile): %s is not a file! Skipping.'
        return 'F'

    # make sure destDir exists
    if not os.path.isdir(destDir):
        try:
            if not dry_run: 
                os.makedirs(destDir)
            if log: 
                print >> log, 'created directory', destDir
        except:
            if log: 
                print >> log, 'failed to create directory', destDir
            return 'F'

    # check if destination file exists (to distinguish between add or
    # update)
    sourceFile = os.path.basename(source)
    destPath = os.path.join(destDir, sourceFile)
    if os.path.isfile(destPath):
        status = 'U'
    else:
        status = 'A'

    # create, log and run the sync command
    syncCmd = 'rsync -uvt'
    if dry_run: syncCmd += 'n'
    cmd = ' '.join([syncCmd, source, destDir])
    if log: 
        print >> log, cmd
    outPipe = os.popen(cmd)
    result = outPipe.readlines()
    if outPipe.close():  # rsync returned non-zero error status
        return 'F'

    # work out return value
    if result[0].strip() != sourceFile:  # no update
        status = 'E'

    return status


### MAIN ###

parser = argparse.ArgumentParser()
parser.add_argument('tmp_dir', help='working directory for downloaded files')
parser.add_argument('target_dir', help='base directory to sort files into')
parser.add_argument('-n', '--dry_run', action="store_true", default=False,
                    help="trial run: write log but don't move files")
args = parser.parse_args()


# download new data into tmp_dir using lftp



# set up for iteration through all files
sourceFiles = []
updatedFiles = []
existingFiles = []
skippedFiles = []
failedFiles = []
LOG = open('sofsUpdate.log', 'w')

print 'sorting files...'
for curDir, dirs, files in os.walk(args.tmp_dir):
    print curDir

    for fileName in files:
        sourcePath = os.path.join(curDir, fileName)
        sourceFiles.append(sourcePath)

        # try to parse filename
        info, err = parseFilename(fileName, minFields=6)

        # if not netcdf file, skip with warning
        if info['extension'] != 'nc':
            print 'WARNING! Non-netCDF file:', fileName
            skippedFiles.append('S '+sourcePath)
            continue

        # work out destination 
        destinationPath = destPath(info, basePath=args.target_dir)
        if not destinationPath:
            print 'WARNING! Unknown destination for', fileName
            skippedFiles.append('S '+sourcePath)
            continue

        # synch file to its destination
        
        status = updateFile(sourcePath, destinationPath, 
                            log=LOG, dry_run=args.dry_run)
        logText = '%s %s' % (status, os.path.join(destinationPath, fileName))
        if status == 'F':
            failedFiles.append(logText)
        elif status == 'E':
            existingFiles.append(logText)
        else:
            updatedFiles.append(logText)


print '%5d files processed' % len(sourceFiles)
print '%5d files updated' % len(updatedFiles)
print '%5d files already at destination' % len(existingFiles)
print '%5d files failed to sync' % len(failedFiles)
print '%5d files skipped' % len(skippedFiles)

print >> LOG, '\nUpdated:\n' + '\n'.join(updatedFiles)
print >> LOG, '\nExisting:\n' + '\n'.join(existingFiles)
print >> LOG, '\nFailed:\n' + '\n'.join(failedFiles)
print >> LOG, '\nSkipped:\n' + '\n'.join(skippedFiles)

LOG.close()
