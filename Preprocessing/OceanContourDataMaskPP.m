function sample_data = imosOceanContourDataMaskPP( sample_data, qcLevel, auto )
%imosOceanContourDataMaskPP apply OceanContour data mask to current variables.
%
% Inputs:
%   sample_data - struct containing the entire data set and dimension data.
%   auto - logical, run QC in batch mode
%
% Outputs:
%   sample_data - same as input, with QC flags added for variable/dimension
%                 data.
%   varChecked  - cell array of variables' name which have been checked
%   paramsLog   - string containing details about params' procedure to include in QC log
%
% Author: Simon Spagnol <s.spagnol@aims.gov.au>
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
narginchk(2, 3);
if ~iscell(sample_data), error('sample_data must be a cell array'); end
if isempty(sample_data), return;                                    end

% no modification of data is performed on the raw FV00 dataset except
% local time to UTC conversion
if strcmpi(qcLevel, 'raw'), return; end


% auto logical in input to enable running under batch processing
if nargin<3, auto=false; end

for k = 1:length(sample_data)

    if isfield(sample_data{k}, 'toolbox_parser') & ~strcmp(sample_data{k}.toolbox_parser, 'OceanContour')
        continue;
    end

    if ~isfield(sample_data{k}.meta, 'twodim_vdatamask')
        continue;
    end
    data_mask = logical(sample_data{k}.meta.twodim_vdatamask);

    varChecked = {};

    % get all necessary dimensions and variables id in sample_data struct
    idUcur = 0;
    idVcur = 0;
    idWcur = 0;
    idWcur2 = 0;
    idCspd = 0;
    idCdir = 0;

    avail_variables = IMOS.get(sample_data{k}.variables, 'name');
    lenVar = length(avail_variables);
    for i=1:lenVar
        paramName = avail_variables{i};
        if strncmpi(paramName, 'UCUR', 4),  idUcur    = i; end
        if strncmpi(paramName, 'VCUR', 4),  idVcur    = i; end
        if strcmpi(paramName, 'WCUR'),      idWcur    = i; end
        if strcmpi(paramName, 'WCUR_2'),    idWcur2   = i; end
        if strcmpi(paramName, 'CSPD'),      idCspd    = i; end
        if strncmpi(paramName, 'CDIR', 4),  idCdir    = i; end
    end

    % check if the data is compatible with the QC algorithm, otherwise quit
    % silently
    idMandatory = (idUcur | idVcur | idWcur | idWcur2 | idCspd | idCdir);
    if ~idMandatory, continue; end

    qcSet = str2double(readProperty('toolbox.qc_set'));
    goodFlag = imosQCFlag('good', qcSet, 'flag');
    badFlag = imosQCFlag('bad', qcSet, 'flag');

    % we try to find out which kind of ADCP we're dealing with and if it is
    % listed as one we should process
    instrument = sample_data{k}.instrument;
    if isfield(sample_data{k}, 'meta')
        if isfield(sample_data{k}.meta, 'instrument_make') && isfield(sample_data{k}.meta, 'instrument_model')
            instrument = [sample_data{k}.meta.instrument_make ' ' sample_data{k}.meta.instrument_model];
        end
    end

    % initially everything is failing the tests
    if idUcur
        idVar = idUcur;
    else
        idVar = idCspd;
    end
    if isfield(sample_data{k}.variables{idVar}, 'flags')
        sizeCur = size(sample_data{k}.variables{idVar}.flags);
    else
        sizeCur = size(sample_data{k}.variables{idVar}.data);
    end
    flags = ones(sizeCur, 'int8')*goodFlag;
    flags(data_mask) = badFlag;

    if idWcur
        sample_data{k}.variables{idWcur}.flags = flags;
        varChecked = [varChecked, {'WCUR'}];
    end

    if idWcur2
        sample_data{k}.variables{idWcur2}.flags = flags;
        varChecked = [varChecked, {'WCUR_2'}];
    end

    if idCspd
        sample_data{k}.variables{idCspd}.flags = flags;
        varChecked = [varChecked, {'CSPD'}];
    end

    if idUcur
        sample_data{k}.variables{idUcur}.flags = flags;
        varChecked = [varChecked, {sample_data{k}.variables{idUcur}.name}];
    end

    if idVcur
        sample_data{k}.variables{idVcur}.flags = flags;
        varChecked = [varChecked, {sample_data{k}.variables{idVcur}.name}];
    end

    if idCdir
        sample_data{k}.variables{idCdir}.flags = flags;
        varChecked = [varChecked, {sample_data{k}.variables{idCdir}.name}];
    end
    
    historyComment = ['OceanContour data mask applied to ' strjoin(varChecked, ',') '.'];

    history = sample_data{k}.history;
    if isempty(history)
        sample_data{k}.history = sprintf('%s - %s', datestr(now_utc, readProperty('exportNetCDF.dateFormat')), ['OceanContourDataMaskPP.m: ' historyComment]);
    else
        sample_data{k}.history = sprintf('%s\n%s - %s', history, datestr(now_utc, readProperty('exportNetCDF.dateFormat')), ['OceanContourDataMaskPP.m: ' historyComment]);
    end
end

end

