function sample_data = magDirVelocityPP( sample_data, qcLevel, auto )
%MAGDIRVELOCITYPP Adds UCUR and VCUR variables to the given data sets, if they
% contain CSPD and CDIR variables.
%
% This function uses trigonometry functions to derive the sea water
% velocity direction and speed data from the meridional and zonal sea water
% velocity speed. It adds the sea water velocity magnitude and direction data 
% as new variables in the data sets. Data sets which do not contain 
% UCUR and VCUR variables or which already contain CSPD and CDIR are left unmodified.
%
% Inputs:
%   sample_data - cell array of data sets, ideally with UCUR and VCUR.
%   qcLevel     - string, 'raw' or 'qc'. Some pp not applied when 'raw'.
%   auto        - logical, run pre-processing in batch mode.
%
% Outputs:
%   sample_data - the same data sets, with CSPD and CDIR variables added.
%
% Based on velocityMagDir.PP
% Author:       Guillaume Galibert <guillaume.galibert@utas.edu.au>
% Author:       Simon Spagnol <s.spagnol@aims.gov.au>

%
% Copyright (C) 2022, Australian Ocean Data Network (AODN) and Integrated 
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
  
  sam = sample_data{k};

  % data set already contains UCUR and VCUR
  if getVar(sam.variables, 'UCUR') && getVar(sam.variables, 'VCUR'), continue; end
  
  cspdIdx = getVar(sam.variables, 'CSPD');
  cdirIdx = getVar(sam.variables, 'CDIR');
  
  % CSPD and CDIR not present in data set
  if ~(cspdIdx && cdirIdx), continue; end
  
  cspd = sam.variables{cspdIdx}.data;
  cdir = sam.variables{cdirIdx}.data;
  
  % convert cspd, cdir CW from N to U/V
  [vcur, ucur] = pol2cart(deg2rad(cdir), cspd);
  
  dimensions = sam.variables{cspdIdx}.dimensions;
  comment = 'magDirVelocityPP.m: UCUR and VCUR were derived from CSPD and CDIR.';
  
  if isfield(sam.variables{cspdIdx}, 'coordinates')
      coordinates = sam.variables{cspdIdx}.coordinates;
  else
      coordinates = '';
  end
    
  % add CSPD and CDIR data as new variable in data set
  sample_data{k} = addVar(...
    sample_data{k}, ...
    'UCUR', ...
    ucur, ...
    dimensions, ...
    comment, ...
    coordinates);

  sample_data{k} = addVar(...
    sample_data{k}, ...
    'VCUR', ...
    vcur, ...
    dimensions, ...
    comment, ...
    coordinates);

    history = sample_data{k}.history;
    if isempty(history)
        sample_data{k}.history = sprintf('%s - %s', datestr(now_utc, readProperty('exportNetCDF.dateFormat')), comment);
    else
        sample_data{k}.history = sprintf('%s\n%s - %s', history, datestr(now_utc, readProperty('exportNetCDF.dateFormat')), comment);
    end
end
