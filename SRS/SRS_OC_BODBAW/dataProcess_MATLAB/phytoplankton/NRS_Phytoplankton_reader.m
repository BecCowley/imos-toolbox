function [DATA,METADATA]=NRS_Phytoplankton_reader(FileToConvert)

% FileToConvert='/home/lbesnard/IMOS/BioOpticalData/TEST_BIO_pigment.csv';
% FileToConvert='/home/lbesnard/IMOS/BioOpticalData/IMOS_FR101997 attribute V2.csv';
% the CSV file is suppose to have three different blocks :'GLOBAL
% ATTRIBUTES', 'TABLE COLUMNS' and 'DATA' .
%
% field delimiter : ,
% Text delimiter : nothing

fid = fopen(FileToConvert);

%% Beginning of 'Global Attributes' Block
tline = fgetl(fid);
while isempty(strfind(tline,'GLOBAL ATTRIBUTES'))
    tline = fgetl(fid);
end

delimiter='|'; %csv delimiter

% nCharacters=length(tline);
% nDelimiter=length(strfind(tline,delimiter));%1 delimiter for 2 fields
% nFields=nDelimiter+1;


%% Metadata
ii=1;
if ~isempty(strfind(tline,'GLOBAL ATTRIBUTES'))
    tline = fgetl(fid);
    
    while isempty(strfind(tline,'TABLE COLUMNS')) && length(strfind(tline,delimiter))~=length(tline) % Condition for the next block
%         disp(tline)
        C = textscan(tline, '%s','Delimiter',delimiter) ;
        AttName{ii}={C{1,1}{1,1}};
        
        % in case there is no attribute name, the cell would have a size
        % equal to 1 only
        if size(C{1},1)>1
            AttVal{ii}={C{1,1}{2,1}};
        else
            AttVal{ii}=[];
        end
        
        ii=ii+1;
        tline = fgetl(fid);
    end
    %     fclose(fid);
else
    disp('bad file format')
    return
end

METADATA=struct;
METADATA.gAttName=AttName;
METADATA.gAttVal=AttVal;

%% Beginning of 'TABLE COLUMNS' Block
% tline = fgetl(fid);
while isempty(strfind(tline,'TABLE COLUMNS'))
    tline = fgetl(fid);
end
%% Variables Attributes
ii=1;

tline = fgetl(fid);
C = textscan(tline, '%s','Delimiter',delimiter);
if strcmpi(C{1,1}{1,1},'Column name') && (strcmpi(C{1,1}{2,1},'CF standard_name') || strcmpi(C{1,1}{2,1},'CF standard name')) && ...
     (strcmpi(C{1,1}{3,1},'IMOS long_name') || strcmpi(C{1,1}{3,1},'IMOS long name')) &&  strcmpi(C{1,1}{4,1},'Units') ...
     && strcmpi(C{1,1}{5,1},'Fill value') && strcmpi(C{1,1}{6,1},'Comments')

    tline = fgetl(fid);
    
    while (isempty(strfind(tline,'DATA')) && length(strfind(tline,delimiter))~=length(tline)) && ~isempty(tline) % Condition for the next block
%           disp(tline)
        C = textscan(tline, '%s','Delimiter',delimiter) ;
        VarName{ii}={C{1,1}{1,1}};
        VarName{ii}=strrep(VarName{ii},'+','_plus_');
        VarName{ii}=strrep(VarName{ii},'-','_');

        
        if size(C{1},1)>1
            VarStandard_Name{ii}={C{1,1}{2,1}};
        else
            VarStandard_Name{ii}=[];
        end
        
        if size(C{1},1)>2
            VarLong_Name{ii}={C{1,1}{3,1}};
        else
            VarLong_Name{ii}=[];
        end
        
        if size(C{1},1)>3
            VarUnits{ii}={C{1,1}{4,1}};
        else
            VarUnits{ii}=[];
        end
        
        if size(C{1},1)>4
            VarFillValue{ii}={C{1,1}{5,1}};
        else
            VarFillValue{ii}=[];
        end
        
        if size(C{1},1)>5
            VarComments{ii}={C{1,1}{6,1}};
        else
            VarComments{ii}=[];
        end
        ii=ii+1;
        tline = fgetl(fid);
    end
else
    disp('bad file format')
    return
end


DATA=struct;
DATA.VarName=VarName;
DATA.Standard_Name=VarStandard_Name;
DATA.Long_Name=VarLong_Name;
DATA.Units=VarUnits;
DATA.FillValue=VarFillValue;
DATA.Comments=VarComments;

% create the Variable line which should be written in the CSV file. This
% line will be check below with tline
VarNameLine=DATA.VarName{1};
for ii=2:length(DATA.VarName)
    VarNameLine=strcat(VarNameLine,',',DATA.VarName{ii});
end

%% Beginning of 'DATA' Block
while isempty(strfind(tline,'DATA'))
    tline = fgetl(fid);
end
%% Var DATA
ii=1;

tline = fgetl(fid);

% if strcmp( tline,VarNameLine)
% tline = fgetl(fid);
C = textscan(tline, '%s','Delimiter',delimiter) ;
VarLineEntries={C{1,1}{:,1}};
VarLineEntries=strrep(VarLineEntries,'+','_plus_');
VarLineEntries=strrep(VarLineEntries,'-','_');
VarLineEntries=VarLineEntries(~cellfun('isempty',VarLineEntries));

if length (VarLineEntries)==length([DATA.VarName{:}])
%     if sum(strcmp(VarLineEntries,[DATA.VarName{:}]))==length(DATA.VarName)
        tline = fgetl(fid);
        
        while ischar(tline) && ~isempty(regexp(tline,'\d', 'once')) && length(strfind(tline,delimiter))~=length(tline) % last block
            %         disp(tline)
            C = textscan(tline, '%s','Delimiter',delimiter) ;
            DataValues{ii}={C{1,1}{:,1}};
            ii=ii+1;
            tline = fgetl(fid);
        end
        %     fclose(fid);
%     else
%         disp('bad file format, the Variable line does not match the with the variable list from the upper block')
%         return
%     end
else
    disp('bad file format, the Variable line does not match the with the variable list from the upper block')
    return
end

fclose(fid);

NVar=length(DATA.VarName);
NlinesData=length(DataValues);
VAL=cell(NVar,NlinesData);
%don't know how to vectorize this
for ii=1:NVar
    for jj=1:NlinesData
        VAL{ii,jj}=DataValues{:,jj}{ii};
    end
end
VAL=VAL';

DATA.Values=VAL;



% {VAL{:,1}}'
% [DATA.VarName{:}]'

clearvars -except DATA METADATA









