function DataFabricFileManagement_MissingQAQC
global FAIMMS_DownloadFolder;
global DataFabricFolder;


LevelName='QAQC';
DataFabricFolder=strcat(DataFabricFolder,'opendap');



if exist(strcat(FAIMMS_DownloadFolder,'/log_archive'),'dir') == 0
    mkdir(strcat(FAIMMS_DownloadFolder,'/log_archive'));
end


%% Copy new files to the DF
LogFilesToCopy=dir(fullfile(FAIMMS_DownloadFolder,strcat('log_ToDo/NoQAQCfile2copy*')));

for tt=1:length(LogFilesToCopy)
    fid2 = fopen(fullfile(FAIMMS_DownloadFolder,'log_ToDo',LogFilesToCopy(tt).name));
    
    kk=1;
    tline = fgetl(fid2);
    while ischar(tline)
        FileTocopy{kk}=tline;
        kk=kk+1;
        tline= fgetl(fid2);
    end
    fclose(fid2);
    
    StatusError=0;
    for kk=1:length(FileTocopy)
        %% create the folder on the DF
        [pathstr] = fileparts(strcat(DataFabricFolder,'/FAIMMS/',FileTocopy{kk}));
        if exist(pathstr,'dir') == 0
            DirCreated=0;
            while ~DirCreated
                DirCreated=mkdir(pathstr);
            end
        end
        [status] = movefile(strcat(FAIMMS_DownloadFolder,filesep,'sorted',filesep,LevelName,filesep,FileTocopy{kk}),strcat(DataFabricFolder,'/FAIMMS/',FileTocopy{kk}));
        if status==0
            StatusError=StatusError+1;
            fprintf('%s - ERROR:  COPY ACHIEVED TO THE DF:  NO --FILE: "%s"\n',datestr(now),strcat(FAIMMS_DownloadFolder,filesep,'sorted',filesep,LevelName,filesep,FileTocopy{kk}));
        elseif status==1
            fprintf('%s - SUCCESS:COPY ACHIEVED TO THE DF: YES --FILE: "%s"\n',datestr(now), strcat(FAIMMS_DownloadFolder,filesep,'sorted',filesep,LevelName,filesep,FileTocopy{kk}));
        end
    end
    
    if StatusError == 0
        movefile(fullfile(FAIMMS_DownloadFolder,'log_ToDo',LogFilesToCopy(tt).name),strcat(FAIMMS_DownloadFolder,'/log_archive'));
    else
        fprintf('%s - ERROR: "%s" has to be check manually,files could not be copied for some reasons\n',datestr(now),LogFilesToCopy(tt).name)
    end
end

%% For any reason if a file has not been copied to the DF, we'll try it again
LevelName='QAQC';


DownloadFolder=strcat(FAIMMS_DownloadFolder,'/sorted/',LevelName);
[~,~,fileNames]=DIRR(DownloadFolder,'.nc','name','isdir','1');


for kk=1:length(fileNames)
    
    %% check it's a file
    if strcmp(fileNames{kk}(end-2:end),'.nc')
        [pathstr, name, ext ~] = fileparts(fileNames{kk});
        
        %% create the folder on the DF
        k=strfind(pathstr,LevelName);
        FileFolder=strcat(DataFabricFolder,'/FAIMMS',pathstr(k+length(LevelName):end));
        if exist(FileFolder,'dir') == 0
            DirCreated=0;%improve dir creation on df : Device or resource busy
            while ~DirCreated
                DirCreated=mkdir(FileFolder);
            end
        end
        
        %% move file to the DF
        [status]  = movefile(fileNames{kk},strcat(FileFolder,filesep,name,ext));
        if status==0
            fprintf('%s - ERROR: "%s" Not copied to the DF:',datestr(now), fileNames{kk});
        end
    end
end

%% Remove old/dupicated files to the DF
LogFilesToDelete=dir(fullfile(FAIMMS_DownloadFolder,strcat('log_ToDo/NoQAQCfile2delete_*')));
for tt=1:length(LogFilesToDelete)
    try
        fid = fopen(fullfile(FAIMMS_DownloadFolder,'log_ToDo',LogFilesToDelete(tt).name));
        
        kk=1;
        tline = fgetl(fid);
        while ischar(tline)
            FileToDelete{kk}=tline;
            kk=kk+1;
            tline= fgetl(fid);
        end
        
        fclose(fid);
        
        SuccessBoolean=1;
        for kk=1:length(FileToDelete)
            if exist(strcat(DataFabricFolder,'/FAIMMS/',FileToDelete{kk}),'file')
                delete(strcat(DataFabricFolder,'/FAIMMS/',FileToDelete{kk}))
                if exist(strcat(DataFabricFolder,'/FAIMMS/',FileToDelete{kk}),'file')
                    fprintf('%s - ERROR: "%s" still exist,can not be deleted\n',datestr(now),strcat(DataFabricFolder,'/FAIMMS/',FileToDelete{kk}))
                    SuccessBoolean=0;
                end
            else
                fprintf('%s - ERROR: "%s" cannot be found for deletion\n',datestr(now),strcat(DataFabricFolder,'/FAIMMS/',FileToDelete{kk}))
                SuccessBoolean=0;
            end
        end
        clear tline
        if SuccessBoolean
            % we move the file from ToDo because everything
            % succeed.Otherwise have to check it manually
            movefile(fullfile(FAIMMS_DownloadFolder,'log_ToDo',LogFilesToDelete(tt).name),strcat(FAIMMS_DownloadFolder,'/log_archive'));
        else
            fprintf('%s - ERROR: "%s" has to be check manually,files could not be deleted for some reasons\n',datestr(now),LogFilesToDelete(tt).name)
        end
        
    catch
        frpintf('%s - No files to delete From the DF',datestr(now))
    end
end