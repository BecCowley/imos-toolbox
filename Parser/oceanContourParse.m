function sample_data = oceanContourParse(filename_in_cell, toolbox_mode)
% function [data] = oceanContourParser(filename_in_cell, toolbox_mode)
%
% The OceanContour Parser for netcdf or mat files.
%
% Inputs:
%
% filename [cell[str]] - A cell containing one filename string.
% toolbox_mode [str] - the processing mode string.
%                    Default: 'timeSeries'
%
% Outputs:
%
% sample_data - toolbox struct containing the data.
%
% Example:
%
% % this is a wrapper,
% % check the `OceanContour.readOceanContourFile`
% % docstring for tests.
%
% author: hugo.oliveira@utas.edu.au
%
narginchk(1, 2)

invalid_file_arg = ~iscellstr(filename_in_cell) && length(filename_in_cell) ~= 1;

if invalid_file_arg
    errormsg('First argument file isn''t a singleton cell of strings')
end

filename = filename_in_cell{1};
inexistent_file = isempty(dir(filename));

if inexistent_file
    errormsg('file %s doesn''t exist', filename)
end

% returns cell array
sample_data = OceanContour.readOceanContourFile(filename);


% check for WAVES files and add any data to the sample_data structure.

[filePath, fileRadName, ~] = fileparts(filename);
fn = strsplit(fileRadName,'.');
txtlist = dir([filePath filesep '*waves.nc']);

if ~isempty(txtlist)
    txtlist = struct2cell(txtlist);
    flist = txtlist(1,:);
    if any(contains(flist, fn{1})) % has a corresponding WAVES file
        waveFile = [filePath, filesep, flist{contains(flist, fn{1})}];
        sample_data = [sample_data, OceanContourWaves.readOceanContourFile(waveFile)];
    end
end
end
