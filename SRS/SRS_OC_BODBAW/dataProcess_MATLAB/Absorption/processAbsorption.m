function processAbsorption(DataFileFolder,xlsFile)
%% processAbsorption
% this processes a XLS absorption file (xlsFile) to CSV and NetCDF files. 
% 2 folders are created, 
% -CSV : the content of this folder has to be copied to the public folder
% of the IMOS cloud storage, at this location SRS/BioOptical
% -NetCDF : the content of this folder has to be copied to the opendap folder
% of the IMOS cloud storage, at this location SRS/BioOptical
% 
% A SQL script is also created to load to the IMOS database in oder to
% populate the table used by geoserver. This script has to be loaded
% manually afterwards.
% Finally, the original XLS file stays at the same location.
%
% Syntax: processAbsorption(DataFileFolder,xlsFile)
%
% Inputs: DATA - structure created by Absorption_CSV_reader
%         METADATA - structure created by Absorption_CSV_reader
%         FileName - filename created by createAbsorptionFilename
%         folderHierarchy - folder structure hierarchy created by createAbsorptionFilename
% Outputs: logfile
%
%
% Example:
%    processAbsorption('/this/is/the/folder','absorptionfile.xls')
%
% Other m-files
% required:
% Other files required:config.txt
% Subfunctions: mkpath
% MAT-files required: none
%
% See also:
% Absorption_CSV_reader,createAbsorptionFilename,processAbsorption,CreateBioOptical_Absorption_SQL_fromCSV,CreateBioOptical_Absorption_NetCDF
%
% Author: Laurent Besnard, IMOS/eMII
% email: laurent.besnard@utas.edu.au
% Website: http://imos.org.au/  http://froggyscripts.blogspot.com
% Aug 2011; Last revision: 28-Nov-2012
Nsheet=str2double(readConfig('xls.absorption.Nsheet', 'config.txt','='));       %number of working sheet in the xls file

CSVfolder=fullfile(DataFileFolder,'CSV');


%% convert xls in csv
% have to convert the xls file manualy in case of bug as a
% csv with OpenOFFICE,UTF8 delimiter column '|' no delimiter text
% each sheet has to be saved in a separate csv file
% need catdoc on linux package sudo apt-get install catdoc
% if there is a bug at this stage, have to save the xls file manualy as a
% csv with OpenOFFICE, delimiter column '|' no delimiter text

% CATDOC version or PERL version. Watch out with 64 bits environment.
% xls2csv_CATDOC_NSheet(filename,Nsheet)
xls2csv_PERL_NSheet([DataFileFolder filesep xlsFile],Nsheet)


%% create the NetCDF file
FilesCSV=dir (strcat(DataFileFolder,filesep,xlsFile(1:end-4),'*.csv'));
for ii=1:length(FilesCSV)
    filename=fullfile(DataFileFolder,char(FilesCSV(ii).name));
    try
        fprintf('%s +++ Process file %s\n',datestr(now),char(FilesCSV(ii).name))
        
        [DATA,METADATA]=Absorption_CSV_reader(filename);
        [FileNameCSV,FileNameNC,folderHierarchy]=createAbsorptionFilename(DATA,METADATA);
        
        % H.2.2. Incomplete multidimensional array representation of time series
        % http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/aphs02.html
        %         [NCfileName,folderHierarchy]=CreateBioOptical_Absorption_NetCDF(DATA,METADATA);
        CreateBioOptical_Absorption_NetCDF2(DATA,METADATA,FileNameNC,folderHierarchy)
        
        %% change csv filename and move to good folder
        mkpath(fullfile(CSVfolder,folderHierarchy))
        movefile(filename,fullfile(CSVfolder,folderHierarchy,char(FileNameCSV)))
        % here we re-write the CSV so the coma is the default
        % delimiter,because some people don't know how do deal with the
        % |'s one .
        csvChangeDelimiter(fullfile(CSVfolder,folderHierarchy,char(FileNameCSV)))
        CreateBioOptical_Absorption_SQL_fromCSV(DATA,METADATA,FileNameCSV,FileNameNC,folderHierarchy)
        plot_absorption_portal2(fullfile(DataFileFolder,filesep,'NetCDF',filesep,folderHierarchy,filesep,char(FileNameNC)),char(FileNameCSV))

    catch
        fprintf('%s - ERROR file %s - REQUIRES DEBUG\n',datestr(now),char(FilesCSV(ii).name))
    end
end


%%% create the PSQL script to load by looking in the NetCDF folder. Not used anymore
% [Files,Bytes,FilesNC] = DIRR(strcat(DataFileFolder,'NetCDF'),'.nc','name');
% for ii=1:length(FilesNC)
%     NCfile=FilesNC{ii};
%     CreateBioOptical_Absorption_SQL(NCfile)
% end
end