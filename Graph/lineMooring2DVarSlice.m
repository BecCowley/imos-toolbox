function lineMooring2DVarSlice(sample_data, varName, yValue, isQC, saveToFile, exportDir)
%LINEMOORING2DVARSLICE Opens a new window from a 2D TimeSeries graph
% displaying variable plotted with  TIME with user selectable value
% of the 2nd dimension.
%
% Inputs:
%   sample_data - struct containing the entire data set and dimension data.
%
%   varName     - string containing the IMOS code for requested parameter.
%
%   timeValue   - double time when the plot must be performed.
%
%   isQC        - logical to plot only good data or not.
%
%   saveToFile  - logical to save the plot on disk or not.
%
%   exportDir   - string containing the destination folder to where the
%               plot is saved on disk.
%
% Author:   Simon Spagnol <s.spagnol@aims.gov.au>
%
% Based on lineMooring2DVarSection by Guillaume Galibert <guillaume.galibert@utas.edu.au>
%

%
% Copyright (C) 2023, Australian Ocean Data Network (AODN) and Integrated
% Marine Observing System (IMOS).
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.
% If not, see <https://www.gnu.org/licenses/gpl-3.0.en.html>.
%
narginchk(6,6);

if ~isstruct(sample_data),  error('sample_data must be a struct');      end
if ~ischar(varName),        error('varName must be a string');          end
if ~isnumeric(yValue),   error('timeValue must be a number');        end
if ~islogical(isQC),        error('isQC must be a logical');            end
if ~islogical(saveToFile),  error('saveToFile must be a logical');      end
if ~ischar(exportDir),      error('exportDir must be a string');        end

monitorRect = getRectMonitor();
iBigMonitor = getBiggestMonitor();

initiateFigure = true;

varTitle = imosParameters(varName, 'long_name');
varUnit  = imosParameters(varName, 'uom');

qcSet = str2double(readProperty('toolbox.qc_set'));

flagList = {'raw', 'good', 'probablyGood', 'probablyBad', 'bad'};
nflags = numel(flagList);
flagVal = zeros([nflags, 1]);
flagDesc = cell(nflags, 1);
flagColor = cell(nflags, 1);
hLineFlag = gobjects(nflags, 1);
% initialize flag values, descriptions, colors
for i = 1:nflags
    name = char(flagList{i});
    value = imosQCFlag(name, qcSet, 'flag');
    flagVal(i) = value;
    flagDesc{i} = imosQCFlag(value, qcSet, 'desc');
    flagColor{i} = imosQCFlag(value, qcSet, 'color');
end

varDesc = cell(3, 1);
hLineVar = gobjects(3, 1);

% instrument description
if ~isempty(strtrim(sample_data.instrument))
    instrumentDesc = sample_data.instrument;
elseif ~isempty(sample_data.toolbox_input_file)
    [~, instrumentDesc] = fileparts(sample_data.toolbox_input_file);
end

if ~isempty(sample_data.meta.depth)
    metaDepth = sample_data.meta.depth;
elseif ~isempty(sample_data.instrument_nominal_depth)
    metaDepth = sample_data.instrument_nominal_depth;
else
    metaDepth = NaN;
end

instrumentSN = '';
if ~isempty(strtrim(sample_data.instrument_serial_number))
    instrumentSN = sample_data.instrument_serial_number;
end

instrumentDesc = [strrep(instrumentDesc, '_', ' ') ' (' num2str(metaDepth) 'm - ' instrumentSN ')'];

%look for TIME and 2nd dimension and relevant variable
iVar = getVar(sample_data.variables, varName);
if iVar == 0
    % something very wrong if we got here
    return;
end
iTime = getVar(sample_data.dimensions, 'TIME');
i2Ddim = sample_data.variables{iVar}.dimensions(2);

[~, indSlice] = min(abs(sample_data.dimensions{i2Ddim}.data - yValue));

dimName = sample_data.dimensions{i2Ddim}.name;
%dimTitle = imosParameters(dimName, 'long_name');
dimUnit  = imosParameters(dimName, 'uom');
dimVals = sample_data.dimensions{i2Ddim}.data;
dimVal = dimVals(indSlice);

xData = datetime(sample_data.dimensions{iTime}.data, 'ConvertFrom', 'datenum');
dataVar = sample_data.variables{iVar}.data;
dataVarSlice = dataVar(:, indSlice);
dimSize = size(sample_data.variables{iVar}.data, 2);

title = [varName ' slice of ' instrumentDesc ' from ' sample_data.deployment_code ' @ ' dimName ' = '  num2str(dimVal) ' ' dimUnit '(index=' num2str(indSlice) ')' ];

backgroundColor = [1 1 1]; % white

if initiateFigure
    fileName = genIMOSFileName(sample_data, '.png');
    visible = 'on';
    if saveToFile, visible = 'off'; end
    visible = 'on';
    if saveToFile, visible = 'off'; end
    hFigVarSlice = figure(...
        'Name',             title, ...
        'NumberTitle',      'off', ...
        'Visible',          visible, ...
        'Color',            backgroundColor, ...
        'OuterPosition',    monitorRect(iBigMonitor, :));
    
    initiateFigure = false;
end

tbh = uitoolbar(hFigVarSlice);
[img, map] = imread(fullfile(matlabroot,'toolbox','matlab','icons','greenarrowicon.gif'));
ptImage = ind2rgb(img, map);
pth_dec = uipushtool(tbh, 'CData', rot90(ptImage, 3), ...
    'Tooltip','Decrease index',...
    'HandleVisibility','on', ...
    'ClickedCallback', @decreaseIndSlice);
pth_inc = uipushtool(tbh, 'CData', rot90(ptImage), ...
    'Tooltip','Increase index', ...
    'HandleVisibility','on', ...
    'ClickedCallback', @increaseIndSlice);

hAxVarSlice = axes;
set(get(hAxVarSlice, 'Title'), 'String', title, 'Interpreter', 'none');
set(get(hAxVarSlice, 'XLabel'), 'String', 'TIME', 'Interpreter', 'none');
set(get(hAxVarSlice, 'YLabel'), 'String', [varTitle ' (' varUnit ')'], 'Interpreter', 'none');

hold(hAxVarSlice, 'on');

% plot data slice
hLineVar(1) = line(xData, dataVarSlice, ...
    'LineStyle', '-');

%varDesc{1} = [varName ' @ ' dimName ' = '  num2str(dimVal) ' ' dimUnit];
varDesc{1} = varName;

% get var QC information
flagsVar = sample_data.variables{iVar}.flags;
flagsVarSlice = flagsVar(:, indSlice);

iBad = flagsVar == flagVal(strcmp('bad', flagList));
iPBad = flagsVar == flagVal(strcmp('probablyBad', flagList));

iFlagsSlice = construct_iFlagSlices(flagsVarSlice, flagVal);

% create nan-ed variable for statistics
iMostlyOk = ~iBad & ~iPBad;
nanData = dataVar;
nanData(~iMostlyOk) = NaN;
[dataVarMean, dataVarStd] = calculate_stats(varName, nanData);

% plot data, mean, std
hLineVar(2) = line(xData, dataVarMean, ...
    'LineStyle',    '-', ... % '--'
    'Color',        [0.9804 0.7843 0.5961]); % pastel orange #FAC898
varDesc{2} = [varName ' mean'];

hLineVar(3) = line(xData, dataVarMean + 3*dataVarStd, ...
    'LineStyle',    '-', ... % '--'
    'Color',        [0.6 0.6 0.6]); % grey

line(xData, dataVarMean - 3*dataVarStd, ...
    'LineStyle',    '-', ...
    'Color',        [0.6 0.6 0.6]); % grey
varDesc{3} = [varName ' mean +/- 3*standard deviation'];

% setup legend entries
for i=1:nflags
    hLineFlag(i) =  line(NaT, NaN, ...
        'LineStyle', 'none', ...
        'Marker', 'o', ...
        'MarkerFaceColor', flagColor{i}, ...
        'MarkerEdgeColor', 'none', ...
        'Visible', 'off'); % this is to make sure all flags are properly displayed within legend
end

% plot flags on top of data profile
for i=1:nflags
    if any(iFlagsSlice(:, i))
        hLineFlag(i) = line(xData(iFlagsSlice(:, i)), dataVarSlice(iFlagsSlice(:, i)), ...
            'LineStyle', 'none', ...
            'Marker', 'o', ...
            'MarkerFaceColor', flagColor{i}, ...
            'MarkerEdgeColor', 'none');
    end
end

% Let's redefine properties after line to make sure grid lines appear
% above color data and XTick and XTickLabel haven't changed
set(hAxVarSlice, ...
    'XGrid',        'on', ...
    'YGrid',        'on', ...
    'Layer',        'top');

% set axes background to be transparent (figure color shows through)
set(hAxVarSlice, 'Color', 'none')

if ~initiateFigure
    iNan = arrayfun(@isempty, hLineVar);
    if any(iNan)
        hLineVar(iNan) = [];
        varDesc(iNan) = [];
    end
    
    hLineVar = [hLineVar; hLineFlag];
    varDesc = [varDesc; flagDesc];
    % Matlab >R2015 legend entries for data which are not plotted
    %   will be shown with reduced opacity.
    % Matlab >R2022 legend will only display plotted variables and
    %  AutoUpdate is ignored.
    legend(hAxVarSlice, ...
        hLineVar,       regexprep(varDesc,'_','\_'), ...
        'Interpreter',  'none', ...
        'Location',     'SouthOutside',...
        'AutoUpdate', 'off');
end

if saveToFile
    fileName = strrep(fileName, '_PLOT-TYPE_', '_LINE-SLICE_'); % IMOS_[sub-facility_code]_[platform_code]_FV01_[time_coverage_start]_[PLOT-TYPE]_C-[creation_date].png
    
    fastSaveas(hFigVarSlice, backgroundColor, fullfile(exportDir, fileName));
    
    close(hFigVarSlice);
end

    function [dataVarMean, dataVarStd] = calculate_stats(varName, nanData)
        % CALCULATE_STATS calculate mean and standard deviation with some
        % handling of directional variables. Note that returned scalar mean
        % and std stats require interpretation.
        
        % Using code derived from circular statistics circstat-matlab package
        % from https://github.com/circstat/circstat-matlab in order to
        % minimize included packages.
        %
        % References:
        %   Statistical analysis of circular data, N. I. Fisher
        %   Topics in circular statistics, S. R. Jammalamadaka et al.
        %   Biostatistical Analysis, J. H. Zar
        %
        %   https://github.com/circstat/circstat-matlab
        
        dirParams = {'CDIR_MAG', 'HEADING_MAG', 'SSDS_MAG', 'SSWD_MAG', ...
            'WPDI_MAG', 'WWPD_MAG', 'SWPD_MAG', ...
            'CDIR', 'HEADING', 'SSDS', 'SSWD', ...
            'WPDI', 'WWPD', 'SWPD'};
        if ismember(varName, dirParams)
            radData = deg2rad(mod(450 - nanData, 360)); % deg CW from N -> rad CCW from E
            
            rsum = sum(exp(complex(0, radData)), 2, 'omitnan');
            
            % mean
            dataVarMean = angle(rsum);
            dataVarMean = mod(450 - rad2deg(dataVarMean), 360); % rad CCW from E -> deg CW from N
            
            % std
            r = abs(rsum)./size(nanData,2);
            dataVarStd = sqrt(2*(1-r)); % angular deviation
            %s0 = sqrt(-2*log(r)); % circular standard deviation
            dataVarStd = mod(450 - rad2deg(dataVarStd), 360);
        else
            dataVarMean = mean(nanData, 2, 'omitnan');
            dataVarStd = std(nanData, 0, 2, 'omitnan');
        end
        
    end

    function iFlagsSlice = construct_iFlagSlices(flagsVarSlice, flag_values)
        
        iFlagsSlice = false([length(flagsVarSlice), numel(flag_values)]);
        for f = 1:numel(flag_values)
            iFlagsSlice(:, f) = flagsVarSlice == flag_values(f);
        end
        % 'probablyBad', 'bad'
        idxPBad = strcmp('probablyBad', flagList);
        idxBad = strcmp('bad', flagList);
        if all(iFlagsSlice(:, idxPBad) | iFlagsSlice(:, idxBad)) && isQC
            [~, fname, fext] = fileparts(sample_data.toolbox_input_file);
            fprintf('%s\n', ['Warning : in ' [fname, fext] ...
                ', there is not any ' varName ' data with good flags at slice ' ...
                ' @ ' dimName ' = '  num2str(dimVal) ' ' dimUnit '(index=' num2str(indSlice) ')']);
        end
        
    end

    function decreaseIndSlice(src, event)
        if indSlice == 1
            return;
        end
        indSlice = max(1, indSlice - 1);
        plotIndSlice(indSlice);
    end

    function increaseIndSlice(src, event)
        if indSlice == dimSize
            return;
        end
        indSlice = min(dimSize, indSlice + 1);
        plotIndSlice(indSlice);
    end

    function plotIndSlice(indSlice)
        h = handle(hLineVar(1));
        h.YData = dataVar(:, indSlice);
        
        dimVal = dimVals(indSlice);
        hAxVarSlice.Title.String = [varName ' slice of ' instrumentDesc ' from ' sample_data.deployment_code ' @ ' dimName ' = '  num2str(dimVal) ' ' dimUnit '(index=' num2str(indSlice) ')' ];
        
        flagsVarSlice = flagsVar(:, indSlice);
        iFlagsSlice = construct_iFlagSlices(flagsVarSlice, flagVal);
        dataVarSlice = dataVar(:, indSlice);
        for k = 1:nflags
            h = handle(hLineFlag(k));
            h.XData = [];
            h.YData = [];
            if any(iFlagsSlice(:, k))
                set(h, 'XData', xData(iFlagsSlice(:, k)), 'YData', dataVarSlice(iFlagsSlice(:, k)));
            end
        end
        refreshdata();
    end

end