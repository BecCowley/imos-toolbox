function [listAllFiles, theoreticalListFiles] = getListFiles(year, month, day, hour, station, theoreticalList)
%This subfunction will return a list of files available on the ARCS
%DATAFABRIC from the date specified and for a specific radar station
%The list of files is stored in the variable "listAllFiles"

%see files radar_WERA_main.m / radar_CODAR_main.m and config.txt for any changes on the
%following global variable
global delayedModeEnd
global delayedMode
global dfradialdata
global dateFormat

listAllFiles = {};
theoreticalListFiles = {};
j = 1;
isComplete = false;

if delayedMode
    [lastYear, lastMonth, lastDay, lastHour, ~, ~] = datevec(delayedModeEnd, dateFormat);
else
    [lastYear, lastMonth, lastDay, lastHour, ~, ~] = datevec(now_utc);
end
lastDate    = datenum(lastYear, lastMonth, lastDay);
lastYear    = num2str(lastYear,     '%i');
lastMonth   = num2str(lastMonth,    '%02i');
lastDay     = num2str(lastDay,      '%02i');

%The search of netCDF files starts at the specified day
listFilesPath = fullfile(dfradialdata, station, year, month, day);
fileInput = fullfile(listFilesPath, '*.nc');
currentDate = now;
listFiles = dir(fileInput);
nFiles = length(listFiles);
for i = 1:nFiles
    % we only want files for that day during and after the specified hour
    % (we actually reprocess the last processed file just in case new radials popped out)
    underScorePos = strfind(listFiles(i).name, '_');
    % we assume the date is the fourth element '_' separated in the file name and has the
    % following format '20120305T042000Z', ex. : IMOS_ACORN_RV_20120305T042000Z_RRK_FV00_radial.nc
    currentHourFile = str2double(listFiles(i).name(underScorePos(3)+10:underScorePos(3)+11));
    if (currentHourFile >= hour)
        curDate = datenum([year, month, day], 'yyyymmdd');
        if curDate == lastDate
            if (currentHourFile <= lastHour)
                % we only consider files that are old enough to be fully copied
                % on disk (older than now - 5min)
                if listFiles(i).datenum + 5/(60*24) < currentDate
                    % we check the file is not corrupted, if it is then is
                    % moved to archive, if not then added to list of files
                    listAllFiles = checkCorrupted(listFilesPath, listFiles(i).name, listAllFiles);
                end
            else
                isComplete = true;
                break;
            end
        else
            % we only consider files that are old enough to be fully copied
            % on disk (older than now - 5min)
            if listFiles(i).datenum + 5/(60*24) < currentDate
                % we check the file is not corrupted, if it is then is
                % moved to archive, if not then added to list of files
                listAllFiles = checkCorrupted(listFilesPath, listFiles(i).name, listAllFiles);
            end
        end
    end
end
clear fileInput listFiles nFiles

if (~isComplete)
    %The search of netCDF files continues at the specified month
    nMonths = fullfile(dfradialdata, station, year, month);
    listDays = dir(nMonths);
    nDays = length(listDays);
    for k = 3:nDays %from 3 to avoid . and ..
        % we only want days after the specified one
        if (str2double(listDays(k).name) > str2double(day))
            curDate = datenum([year, month, listDays(k).name], 'yyyymmdd');
            if curDate <= lastDate
                if curDate == lastDate
                    % we take all the files for each of those days until the
                    % last day, last hour
                    listFilesPath = fullfile(nMonths, listDays(k).name);
                    fileInput = fullfile(listFilesPath, '*.nc');
                    listFiles = dir(fileInput);
                    nFiles = length(listFiles);
                    for i = 1:nFiles
                        underScorePos = strfind(listFiles(i).name, '_');
                        currentHourFile = str2double(listFiles(i).name(underScorePos(3)+10:underScorePos(3)+11));
                        if (currentHourFile <= lastHour)
                            % we only consider files that are old enough to be fully copied
                            % on disk (older than now - 5min)
                            if listFiles(i).datenum + 5/(60*24) < currentDate
                                % we check the file is not corrupted, if it is then is
                                % moved to archive, if not then added to list of files
                                listAllFiles = checkCorrupted(listFilesPath, listFiles(i).name, listAllFiles);
                            end
                        else
                            isComplete = true;
                            break;
                        end
                    end
                    if isComplete, break; end
                else
                    % we take all the files for each of those days
                    listFilesPath = fullfile(nMonths, listDays(k).name);
                    fileInput = fullfile(listFilesPath, '*.nc');
                    listFiles = dir(fileInput);
                    nFiles = length(listFiles);
                    for i = 1:nFiles
                        % we only consider files that are old enough to be fully copied
                        % on disk (older than now - 5min)
                        if listFiles(i).datenum + 5/(60*24) < currentDate
                            % we check the file is not corrupted, if it is then is
                            % moved to archive, if not then added to list of files
                            listAllFiles = checkCorrupted(listFilesPath, listFiles(i).name, listAllFiles);
                        end
                    end
                end
            else
                isComplete = true;
                break;
            end
        end
    end
    clear nMonths listDays nDays fileInput listFiles nFiles
end

if (~isComplete)
    %The search of netCDF files continues at the specified year
    yearInput = fullfile(dfradialdata, station, year);
    listMonths = dir(yearInput);
    nMonth = length(listMonths);
    for k = 3:nMonth %from 3 to avoid . and ..
        % we only want months after the specified one
        if (str2double(listMonths(k).name) > str2double(month))
            % we take all the days for each of those months
            nMonths = fullfile(yearInput, listMonths(k).name);
            listDays = dir(nMonths);
            nDays = length(listDays);
            for l = 3:nDays
                curDate = datenum([year, listMonths(k).name, listDays(l).name], 'yyyymmdd');
                if curDate <= lastDate
                    if curDate == lastDate
                        % we take all the files for each of those days until the
                        % last day, last hour
                        listFilesPath = fullfile(nMonths, listDays(l).name);
                        fileInput = fullfile(listFilesPath, '*.nc');
                        listFiles = dir(fileInput);
                        nFiles = length(listFiles);
                        for i = 1:nFiles
                            underScorePos = strfind(listFiles(i).name, '_');
                            currentHourFile = str2double(listFiles(i).name(underScorePos(3)+10:underScorePos(3)+11));
                            if (currentHourFile <= lastHour)
                                % we only consider files that are old enough to be fully copied
                                % on disk (older than now - 5min)
                                if listFiles(i).datenum + 5/(60*24) < currentDate
                                    % we check the file is not corrupted, if it is then is
                                    % moved to archive, if not then added to list of files
                                    listAllFiles = checkCorrupted(listFilesPath, listFiles(i).name, listAllFiles);
                                end
                            else
                                isComplete = true;
                                break;
                            end
                        end
                        if isComplete, break; end
                    else
                        % we take all the files for each of those days
                        listFilesPath = fullfile(nMonths, listDays(l).name);
                        fileInput = fullfile(listFilesPath, '*.nc');
                        listFiles = dir(fileInput);
                        nFiles = length(listFiles);
                        for i = 1:nFiles
                            % we only consider files that are old enough to be fully copied
                            % on disk (older than now - 5min)
                            if listFiles(i).datenum + 5/(60*24) < currentDate
                                % we check the file is not corrupted, if it is then is
                                % moved to archive, if not then added to list of files
                                listAllFiles = checkCorrupted(listFilesPath, listFiles(i).name, listAllFiles);
                            end
                        end
                    end
                else
                    isComplete = true;
                    break;
                end
                if isComplete, break; end
            end
        end
    end
    clear yearInput listMonths nMonth nMonths listDays nDays fileInput listFiles nFiles
end

if (~isComplete)
    %Let's finally look for netCDF files for the following years after the
    %specified one
    stationInput = fullfile(dfradialdata, station);
    listYears = dir(stationInput);
    nYears = length(listYears);
    for k = 3:nYears %from 3 to avoid . and ..
        % we only want years after the specified one
        if (str2double(listYears(k).name) > str2double(year))
            % we take all the months for each of those years
            yearInput = fullfile(stationInput, listYears(k).name);
            listMonths = dir(yearInput);
            nMonth = length(listMonths);
            for m = 3:nMonth
                % we take all the days for each of those months
                monthInput = fullfile(yearInput, listMonths(m).name);
                listDays = dir(monthInput);
                nDays = length(listDays);
                for l = 3:nDays
                    curDate = datenum([listYears(k).name, listMonths(m).name, listDays(l).name], 'yyyymmdd');
                    if curDate <= lastDate
                        if curDate == lastDate
                            % we take all the files for each of those days until the
                            % last day, last hour
                            listFilesPath = fullfile(monthInput, listDays(l).name);
                            fileInput = fullfile(listFilesPath, '*.nc');
                            listFiles = dir(fileInput);
                            nFiles = length(listFiles);
                            for i = 1:nFiles
                                underScorePos = strfind(listFiles(i).name, '_');
                                currentHourFile = str2double(listFiles(i).name(underScorePos(3)+10:underScorePos(3)+11));
                                if (currentHourFile <= lastHour)
                                    % we only consider files that are old enough to be fully copied
                                    % on disk (older than now - 5min)
                                    if listFiles(i).datenum + 5/(60*24) < currentDate
                                        % we check the file is not corrupted, if it is then is
                                        % moved to archive, if not then added to list of files
                                        listAllFiles = checkCorrupted(listFilesPath, listFiles(i).name, listAllFiles);
                                    end
                                else
                                    isComplete = true;
                                    break;
                                end
                            end
                            if isComplete, break; end
                        else
                            % we take all the files for each of those days
                            listFilesPath = fullfile(monthInput, listDays(l).name);
                            fileInput = fullfile(listFilesPath, '*.nc');
                            listFiles = dir(fileInput);
                            nFiles = length(listFiles);
                            for i = 1:nFiles
                                % we only consider files that are old enough to be fully copied
                                % on disk (older than now - 5min)
                                if listFiles(i).datenum + 5/(60*24) < currentDate
                                    % we check the file is not corrupted, if it is then is
                                    % moved to archive, if not then added to list of files
                                    listAllFiles = checkCorrupted(listFilesPath, listFiles(i).name, listAllFiles);
                                end
                            end
                        end
                    else
                        isComplete = true;
                        break;
                    end
                    if isComplete, break; end
                end
                if isComplete, break; end
            end
        end
    end
    clear stationInput listYears nYears yearInput listMonths nMonth monthInput listDays nDays fileInput listFiles nFiles
end

if ~isempty(listAllFiles)
    switch station
        case {'LEI', 'CSP', 'FRE', 'NNB'}
            minute = 5;
        otherwise
            minute = 0;
    end
    
    % check that the station's minute hasn't been swapped from 0 to 5 (or
    % the other)
    underScorePos = cell2mat(strfind(listAllFiles, '_'));
    charListAllFiles = char(listAllFiles);
    dateListAllFiles = charListAllFiles(:, underScorePos(:,3)+1:underScorePos(:,3)+16);
    yearListAllFiles = str2num(dateListAllFiles(:, 1:4)); % str2num is used here because we deal with chars
    monthListAllFiles = str2num(dateListAllFiles(:, 5:6));
    dayListAllFiles = str2num(dateListAllFiles(:, 7:8));
    minuteListAllFiles = str2num(dateListAllFiles(:, 13));
    
    iSwappedMinute = (minuteListAllFiles == 5 - minute);
    if any(iSwappedMinute) && theoreticalList
        impactedYear = yearListAllFiles(iSwappedMinute);
        impactedMonth = monthListAllFiles(iSwappedMinute);
        impactedDay = dayListAllFiles(iSwappedMinute);
        impactedDate = datenum(impactedYear, impactedMonth, impactedDay);
        impactedDate = unique(impactedDate);
        for i=1:length(impactedDate)
            disp(['Warning : ' station ' files have their minute swapped (not consistent with the minute ' num2str(minute) ') for the date ' datestr(impactedDate(i), 'yyyy-mm-dd')]);
        end
    end

    % what is best for the following processing is to have a theoretical list
    % of files if we would have all the files found on the DF
    startDate   = datestr(datenum(str2double(year), str2double(month), str2double(day), hour, minute, 0), dateFormat);
    
    lastFile = listAllFiles{end};
    
    underScorePos = strfind(lastFile, '_');
    % we assume the date is the fourth element '_' separated in the file name and has the
    % following format '20120305T042000Z', ex. : IMOS_ACORN_RV_20120305T042000Z_RRK_FV00_radial.nc
    endDate = lastFile(underScorePos(3)+1:underScorePos(3)+16);
    endDate(12) = '5';
    endDate(13) = num2str(minute);
    
    nIdealFiles = round((datenum(endDate, dateFormat) - datenum(startDate, dateFormat))/((1/24)/6)) + 1;
    allDatesNum = (datenum(startDate, dateFormat):((1/24)/6):datenum(endDate, dateFormat))';
    
    prefix = repmat(lastFile(1:underScorePos(3)), nIdealFiles, 1);
    suffix = repmat(lastFile(underScorePos(4)-1:end), nIdealFiles, 1);
    
    theoreticalListFiles = cellstr([prefix datestr(allDatesNum, dateFormat) suffix]);
    
    if theoreticalList
        % we set to empty where there is no file
        iExist = ismember(theoreticalListFiles, listAllFiles);
        listAllFiles = theoreticalListFiles;
        if any(~iExist), [listAllFiles{~iExist}] = deal(''); end
    end
end
end

function listAllFiles = checkCorrupted(listFilesPath, ncFile, listAllFiles)
% we check ncFile is not corrupted, if it is then is
% moved to archive, if not then added to listAllFiles
                                    
currentFilePath = fullfile(listFilesPath, ncFile);
status = system(['ncdump ' currentFilePath ' &> /dev/null']);
if (status == 0)
    listAllFiles{end+1, 1} = ncFile;
else
    status = system(['rsync -aq --remove-source-files ' currentFilePath ' ' getenv('ARCHIVE') '/ACORN/corrupted/']);
    if (status == 0)
        fprintf('%s\r\n', ['Corrupted file ' currentFilePath ' moved to $ARCHIVE/ACORN/corrupted/']);
    else
        fprintf('%s\r\n', ['Corrupted file ' currentFilePath ' could not be moved to $ARCHIVE/ACORN/corrupted/']);
    end
end
end