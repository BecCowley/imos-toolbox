classdef OceanContourWaves
    %classdef OceanContour
    %
    % This is a class containing methods that defines several
    % fields and functions related to the OceanContour Parser.
    % This includes utility functions and variable/attribute
    % mappings to the toolbox structures.
    %
    % author: hugo.oliveira@utas.edu.au
    %
    %TODO: Design metadata typecasting.
    properties (Constant)
        beam_angles = struct('Signature55', 20, 'Signature100', 20, 'Signature250', 20, 'Signature500', 25, 'Signature1000', 25);
    end
    
    methods (Static)
        
        function metaname = build_meta_attr_midname(group_name)
            %function metaname = build_meta_attr_midname(group_name)
            %
            % Generate the middle name for global attributes given
            % a group name.
            %
            % Input:
            %
            % group_name - the fieldname or netcdf group name.
            %
            % Output:
            %
            % metaname - the mid/partial metadata attribute name string.
            %
            % Example:
            %
            % midname = OceanContourWaves.build_meta_attr_midname('Avg');
            % assert(strcmp(midname,'avg'))
            % midname = OceanContourWaves.build_meta_attr_midname('burstAltimeter');
            % assert(strcmp(midname,'burstAltimeter'))
            %
            if ~ischar(group_name)
                errormsg('first argument is not a string')
            end
            
            metaname = [lower(group_name(1)) group_name(2:end)];
        end
        
        function attname = build_instrument_name(group_name, var_name)
            %function attname = build_instrument_name(group_name,var_name)
            %
            % Generate instrument tokens for the attribute names
            % for in OceanContour files.
            %
            % The token is a three part string:
            % part1 - "Instrument" string, followed
            % part2 -  group/dataset name (with the first  letter lower)
            % part3 - the "variable" token/name.
            %
            % Inputs:
            %
            % group_name [str] - the dataset group (field) name
            % var_name [str] - the variable name (last token).
            %
            % Output:
            % attname [str] - the attribute name.
            %
            % Example:
            % name = OceanContourWaves.build_instrument_name('Avg', 'coordSystem');
            % assert(strcmpi(name,'Instrument_avg_coordSystem'))
            %
            %
            % author: hugo.oliveira@utas.edu.au
            %
            narginchk(2, 2)
            
            if ~ischar(group_name)
                errormsg('first argument is not a string')
            elseif ~ischar(var_name)
                errormsg('second argument is not a string')
            end
            
            meta_attr_midname = OceanContourWaves.build_meta_attr_midname(group_name);
            attname = ['Instrument_' meta_attr_midname '_' var_name];
        end
        
        
        function verify_mat_groups(matdata)
            %just raise a proper error for invalid OceanContour mat files.
            try
                matdata.Config;
            catch
                errormsg('%s do not contains the ''Config'' metadata fieldname', filename)
            end
            
            ngroups = numel(fieldnames(matdata));
            
            if ngroups < 2
                errormsg('%s do not contains any data fieldname', fielname)
            end
            
        end
        
        function verify_netcdf_groups(info)
            %just raise a proper error for invalid OceanContour netcdf groups.
            try
                assert(strcmp(info.Groups(1).Name, 'Config'))
                assert(strcmp(info.Groups(2).Name, 'Data'))
            catch
                errormsg('contains an invalid OceanContour structure. please report this error with your data file: %s', filename)
            end
            
        end
        
        function warning_failed(failed_items, filename)
            %just raise a proper warning for failed variable reads
            for k = 1:numel(failed_items)
                dispmsg('%s: Couldn''t read variable `%s` in %s', mfilename, failed_items{k}, filename)
            end
            
        end
        
        function [attmap] = get_attmap(file_metadata, ftype, group_name)
            %function [attmap] = get_attmap(ftype, group_name)
            %
            % Generate dynamical attribute mappings based on
            % the dataset group name.
            %
            % Inputs:
            %
            % file_metadata []
            % ftype [str] - the file type. 'mat' or 'netcdf';
            % group_name [str] - the OceanContour dataset group name.
            %
            % Outputs:
            %
            % attmap [struct[str,str]] - mapping between imos attributes
            %                           & OceanContour attributes.
            %
            %
            % Example:
            %
            % %basic usage
            % attmap = OceanContourWaves.get_attmap('Avg');
            % fnames = fieldnames(attmap);
            % assert(contains(fnames,'instrument_model'))
            % original_name =attmap.instrument_model;
            % assert(strcmp(original_name,'Instrument_instrumentName'));
            %
            % author: hugo.oliveira@utas.edu.au
            %
            
            if ~ischar(ftype)
                errormsg('First argument is not a string')
            elseif ~strcmpi(ftype, 'mat') && ~strcmpi(ftype, 'netcdf')
                errormsg('First argument %s is an invalid ftype. Accepted file types are ''mat'' and ''netcdf''.', ftype)
            elseif ~ischar(group_name)
                errormsg('Second argument is not a string')
            end
            
            attmap = struct();
            for k = fieldnames(file_metadata)'
                key = k{1};
                attmap.(key) = file_metadata.(key);
            end
            meta_attr_midname = OceanContourWaves.build_meta_attr_midname(group_name);
            
            attmap.('instrument_model') = 'Instrument_instrumentName';
            attmap.('beam_angle') = 'DataInfo_slantAngles';
            %attmap.('beam_interval') = 'DataInfo_slantAngles'; % what is this?
            attmap.('coordinate_system') = OceanContourWaves.build_instrument_name(group_name, 'coordSystem');
            attmap.('converted_to_enu') = 'DataInfo_transformsAndCorrections_addENU';
            
            % while this might be a waves file, some info is 'avg' prefix
            avg_group_name = 'avg';
            attmap.('nBeams') = OceanContourWaves.build_instrument_name(avg_group_name, 'nBeams');
            attmap.('activeBeams') = OceanContourWaves.build_instrument_name(avg_group_name, 'activeBeams'); %no previous name
            attmap.('magDec_User') = 'Instrument_user_decl';
            attmap.('magDec_DataInfo') = 'DataInfo_transformsAndCorrections_magneticDeclination';
            attmap.('binMapping') = 'DataInfo_transformsAndCorrections_binMapping';
            attmap.('binMapping_applied') = 'DataInfo_transformsAndCorrections_binMapping_description';
            
            if strcmpi(ftype, 'mat') || strcmpi(ftype, 'netcdf')
                try
                    attmap.('instrument_serial_no') = 'Instrument_serialNumberDoppler';
                    attmap.('binSize') = OceanContourWaves.build_instrument_name(avg_group_name, 'cellSize');
                catch
                end
            end
            
            %custom & dynamical fields
            attmap.(['instrument_' meta_attr_midname '_enable']) = OceanContourWaves.build_instrument_name(avg_group_name, 'enable');
            
            switch meta_attr_midname
                case 'avg'
                    attmap.('instrument_avg_interval') = OceanContourWaves.build_instrument_name(group_name, 'averagingInterval');
                    attmap.('instrument_sample_interval') = OceanContourWaves.build_instrument_name(group_name, 'measurementInterval');
                    %TODO: need a more complete file to test below below
                case 'burst'
                    attmap.('instrument_burst_interval') = OceanContourWaves.build_instrument_name(group_name, 'burstInterval');
                case 'bursthr'
                    attmap.('instrument_bursthr_interval') = OceanContourWaves.build_instrument_name(group_name, 'burstHourlyInterval');
                case 'burstAltimeter'
                    attmap.('instrument_burstAltimeter_interval') = OceanContourWaves.build_instrument_name(group_name, 'burstAltimeterInterval');
                case 'burstRawAltimeter'
                    attmap.('instrument_burstRawAltimeter_interval') = OceanContourWaves.build_instrument_name(group_name, 'burstRawAltimeterInterval');
                case 'waves'
                    attmap.('instrument_sample_interval') = OceanContourWaves.build_instrument_name('burst', 'measurementInterval');
                    attmap.('nSamples') = OceanContourWaves.build_instrument_name('burst', 'nSamples');
                    attmap.('sampleRate') = OceanContourWaves.build_instrument_name('burst', 'sampleRate');
                    attmap.('coordinate_system') = OceanContourWaves.build_instrument_name('avg', 'coordSystem');
            end
            
        end
        
        function [varmap] = get_varmap(ftype, group_name, nbeams, custom_magnetic_declination, binmapped)
            %function [varmap] = get_varmap(ftype, group_name,nbeams,custom_magnetic_declination)
            %
            % Generate dynamical variable mappings for a certain
            % group of variables, given the number of beams and if custom
            % magnetic adjustments were made.
            %
            % Inputs:
            %
            % ftype [str] - The file type. 'mat' or 'netcdf'.
            % group_name [str] - the OceanContour dataset group name.
            % nbeams [double] - The nbeams used on the dataset.
            % custom_magnetic_declination [logical] - true for custom
            %                                         magnetic values.
            %
            % Outputs:
            %
            % vttmap [struct[str,str]] - mapping between imos variables
            %                           & OceanContour variables.
            %
            %
            % Example:
            %
            % %basic usage
            %
            % varmap = OceanContourWaves.get_attmap('Avg',4,False);
            % assert(strcmp(attmap.WCUR_2,'Vel_Up2'));
            %
            % % nbeams == 3
            % varmap = OceanContourWaves.get_varmap('Avg',3,False);
            % f=false;try;varmap.WCUR_2;catch;f=true;end
            % assert(f)
            %
            % % custom magdec - may change with further testing
            % varmap = OceanContourWaves.get_varmap('Avg',4,True);
            % assert(strcmp(varmap.UCUR_MAG,'Vel_East'))
            %
            %
            % author: hugo.oliveira@utas.edu.au
            %
            narginchk(5, 5)
            
            if ~ischar(ftype)
                errormsg('First argument is not a string')
            elseif ~strcmpi(ftype, 'mat') && ~strcmpi(ftype, 'netcdf')
                errormsg('First argument %s is an invalid ftype. Accepted file types are ''mat'' and ''netcdf''.', ftype)
            elseif ~ischar(group_name)
                errormsg('Second argument is not a string')
            elseif ~isscalar(nbeams)
                errormsg('Third argument is not a scalar')
            elseif ~islogical(custom_magnetic_declination)
                errormsg('Fourth argument is not logical')
            elseif ~islogical(binmapped)
                errormsg('Fifth argument is not logical')
            end
            
            is_netcdf = strcmpi(ftype, 'netcdf');
            
            json_varwaves = jsondecode(fileread('oceancontour_waves_variables.json'));
            json_varmap = jsondecode(fileread('oceancontour_waves_map.json'));
            
            varmap = struct();
            varmap.('binSize') = 'CellSize';
            varmap.('TIME') = 'MatlabTimeStamp';
            
            prefix = '';
            if binmapped
                prefix = 'BinMap';
            end
            
            if is_netcdf
                varmap.('instrument_serial_no') = 'SerialNumber';
                %TODO: reinforce uppercase at first letter? nEed to see more files.
                varmap.('HEIGHT_ABOVE_SENSOR') = [prefix group_name 'VelocityENU_Range'];
                %TODO: Handle magnetic & along beam cases.
                %varmap.('DIST_ALONG_BEAMS') = [group_name 'Velocity???_Range'];
                %TODO: evaluate if when magnetic declination is provided, the
                %velocity fields will be corrected or not (as well as any rename/comments added).
                varmap.(ucur_name) = [prefix 'Vel_East'];
                varmap.(vcur_name) = [prefix 'Vel_North'];
                varmap.(heading_name) = 'Heading';
                varmap.('WCUR') = [prefix 'Vel_Up1'];
                varmap.('ABSI1') = [prefix 'Amp_Beam1'];
                varmap.('ABSI2') = [prefix 'Amp_Beam2'];
                varmap.('ABSI3') = [prefix 'Amp_Beam3'];
                varmap.('CMAG1') = [prefix 'Cor_Beam1'];
                varmap.('CMAG2') = [prefix 'Cor_Beam2'];
                varmap.('CMAG3') = [prefix 'Cor_Beam3'];
                
                if nbeams > 3
                    varmap.('WCUR_2') = [prefix 'Vel_Up2'];
                    varmap.('ABSI4') = [prefix 'Amp_Beam4'];
                    varmap.('CMAG4') = [prefix 'Cor_Beam4'];
                end
                
            else
                %TODO: check if norteks also change the variable names
                %when exporting to matlab.
                %instrument_serial_no is on metadata for matfiles.
                varmap.('HEIGHT_ABOVE_SENSOR') = 'Range';
                varmap.(ucur_name) = 'VelEast';
                varmap.(vcur_name) = 'VelNorth';
                varmap.(heading_name) = 'Heading';
                varmap.('WCUR') = 'VelUp1';
                varmap.('ABSI1') = 'AmpBeam1';
                varmap.('ABSI2') = 'AmpBeam2';
                varmap.('ABSI3') = 'AmpBeam3';
                varmap.('CMAG1') = 'CorBeam1';
                varmap.('CMAG2') = 'CorBeam2';
                varmap.('CMAG3') = 'CorBeam3';
                
                if nbeams > 3
                    varmap.('WCUR_2') = 'VelUp2';
                    varmap.('ABSI4') = 'AmpBeam4';
                    varmap.('CMAG4') = 'CorBeam4';
                end
                
            end
            
            varmap.('data_mask') = 'DataMask';
            varmap.('status') = 'Status';
            varmap.('TEMP') = 'WaterTemperature';
            varmap.('PRES_REL') = 'Pressure';
            varmap.('SSPD') = 'SpeedOfSound';
            varmap.('BAT_VOLT') = 'Battery';
            varmap.('PITCH') = 'Pitch';
            varmap.('ROLL') = 'Roll';
            varmap.('ERROR') = 'Error';
            varmap.('AMBIG_VEL') = 'Ambiguity';
            varmap.('TRANSMIT_E') = 'TransmitEnergy';
            varmap.('NOMINAL_CORR') = 'NominalCor';
            
            
        end
        
       
        function imos_name = get_imos_mapped_name(var_name, var_map, mag_dec)
            % TODO: if a directional variable has not undergone
            % magnetic declination append '_MAG', but is it the variable
            % and or the dimension that should be fixed?
            mag_params = {'CurrentDirection' 'Direction_DirTp'...
                'Direction_SprTp' 'Direction_MeanDir' 'Heading'...
                'Direction' 'ASTSpectra_Direction' 'PressureSpectra_Direction'...
                'VelocitySpectra_Direction'};
            imos_name = var_name;
            mapped_name = var_map(var_name).imos_name;
            if ~isempty(mapped_name)
                imos_name = mapped_name;
            end
            % Not implemented yet until magneticDeclinationPP can be
            % updated to handled required transforms
%             if ~mag_dec
%                 if ismember(var_name, mag_params)
%                     imos_name = strcat(imos_name, '_MAG');
%                 end
%             end
                
        end
        
        
        function [sample_data] = readOceanContourFile(filename)
            % function [sample_data] = readOceanContourFile(filename)
            %
            % Read an OceanContour netcdf or mat file and convert fields
            % to the matlab toolbox structure. Variables are read
            % as is.
            %
            % Supported Innstruments: Nortek ADCP Signatures.
            %
            % The Ocean contour software write nested netcdf4 groups:
            % > root
            %    |
            % {root_groups}
            %    | -> Config ["global" file metadata only]
            %    | -> Data [file datasets leaf]
            %          |
            %       {data_groups}
            %          | -> Avg [data+metadata]
            %          | -> ... [data+metadata]
            %
            % Or a flat mat file:
            % > root
            %    |
            %{data_groups}
            %      | -> [dataset-name] [data]
            %      | -> Config [metadata]
            %
            %
            % Inputs:
            %
            % filename [str] - the filename.
            %
            % Outputs:
            %
            % sample_data - the toolbox structure.
            %
            % Example:
            %
            % %read from netcdf
            % file = [toolboxRootPath 'data/testfiles/netcdf/Nortek/OceanContour/Signature/s500_enu_avg.nc'];
            % [sample_data] = OceanContourWaves.readOceanContourFile(file);
            % assert(strcmpi(sample_data{1}.meta.instrument_model,'Signature500'))
            % assert(isequal(sample_data{1}.meta.instrument_avg_interval,60))
            % assert(isequal(sample_data{1}.meta.instrument_sample_interval,600))
            % assert(strcmpi(sample_data{1}.meta.coordinate_system,'ENU'))
            % assert(isequal(sample_data{1}.meta.nBeams,4))
            % assert(strcmpi(sample_data{1}.dimensions{2}.name,'HEIGHT_ABOVE_SENSOR'))
            % assert(~isempty(sample_data{1}.variables{end}.data))
            %
            % % read from matfile
            % file = [toolboxRootPath 'data/testfiles/mat/Nortek/OceanContour/Signature/s500_enu_avg.mat'];
            % [sample_data] = OceanContourWaves.readOceanContourFile(file);
            % assert(strcmpi(sample_data{1}.meta.instrument_model,'Signature500'))
            % assert(isequal(sample_data{1}.meta.instrument_avg_interval,60))
            % assert(isequal(sample_data{1}.meta.instrument_sample_interval,600))
            % assert(strcmpi(sample_data{1}.meta.coordinate_system,'ENU'))
            % assert(isequal(sample_data{1}.meta.nBeams,4))
            % assert(strcmpi(sample_data{1}.dimensions{2}.name,'HEIGHT_ABOVE_SENSOR'))
            % assert(~isempty(sample_data{1}.variables{end}.data))
            %
            %
            % author: hugo.oliveira@utas.edu.au
            %
            narginchk(1, 1)
            
            try
                info = ncinfo(filename);
                ftype = 'netcdf';
                
            catch
                % TODO: handle matlab wave exported file
%                 errormsg('%s is not netcdf file', filename)
%                 try
%                     matdata = load(filename);
%                     ftype = 'mat';
%                 catch
%                     errormsg('%s is not a mat or netcdf file', filename)
%                 end
                
            end
            
            is_netcdf = strcmpi(ftype, 'netcdf');
            
            if is_netcdf
                OceanContourWaves.verify_netcdf_groups(info);
                file_metadata = nc_flat(info.Groups(1).Attributes, false);
                data_metadata = nc_flat(info.Groups(2).Groups, false);
                
                ncid = netcdf.open(filename);
                root_groups = netcdf.inqGrps(ncid);
                data_group = root_groups(2);
                
                dataset_groups = netcdf.inqGrps(data_group);
                get_group_name = @(x)(netcdf.inqGrpName(x));
                
            else
                OceanContourWaves.verify_mat_groups(matdata);
                file_metadata = matdata.Config;
                matdata = rmfield(matdata, 'Config'); %mem optimisation.
                
                dataset_groups = fieldnames(matdata);
                get_group_name = @(x)(getindex(split(x, '_Data'), 1));
                
            end
            
            n_datasets = numel(dataset_groups);
            sample_data = cell(1, n_datasets);
            
            json_varwaves = jsondecode(fileread('oceancontour_waves_variables.json'));
            json_varmap = jsondecode(fileread('oceancontour_waves_map.json'));
            var_map = containers.Map;
            var_mapping = struct();
            for i = 1:numel(json_varmap.mapping)
                vname = json_varmap.mapping(i).name;
                var_map(vname) = json_varmap.mapping(i);
                var_mapping.(vname) = json_varmap.mapping(i).name;
            end
            
            known_var_names = {json_varwaves.variables.name};
            wave_dim_names = {json_varwaves.dimensions.name};
                
            % TODO: merge into OceanContour.m to parse either avg or burst
            % file. At the moment this is waves only.
            for k = 1:n_datasets
                
                % start by loading preliminary information into the metadata struct, so we
                % can define the variable names and variables to import.
                meta = struct();
                
                group_name = get_group_name(dataset_groups);
                meta_attr_midname = OceanContourWaves.build_meta_attr_midname(group_name);
                
                %load toolbox_attr_names:file_attr_names dict.
                att_mapping = OceanContourWaves.get_attmap(file_metadata, ftype, group_name);
                
                %access pattern - use lookup based on expected names,
                get_att = @(x)(file_metadata.(att_mapping.(x)));
                
                nBeams = double(get_att('nBeams'));
                
                try
                    activeBeams = double(get_att('activeBeams'));
                catch
                    activeBeams = Inf;
                end
                
                meta.nBeams = min(nBeams, activeBeams);
                
                %                 try
                %                     assert(meta.nBeams == 4);
                %                     %TODO: support variable nBeams. need more files.
                %                 catch
                %                     errormsg('Only 4 Beam ADCP are supported. %s got %d nBeams', filename, meta.nBeams)
                %                 end
                
                magDec_User = get_att('magDec_User');
                magDec_DataInfo = get_att('magDec_DataInfo');
                has_magdec_user = logical(magDec_User);
                has_magdec_oceancontour = logical(magDec_DataInfo);
                meta.magDec = 0.0;
                custom_magnetic_declination = has_magdec_user | has_magdec_oceancontour;
                if has_magdec_oceancontour
                    meta.magDec = magDec_DataInfo;
                elseif has_magdec_user
                    meta.magDec = magDec_User;
                end
                
                binmapped = false;
                try
                    meta.binMapping = get_att('binMapping');
                    binmapped = logical(meta.binMapping);
                catch
                    binmapped = false;
                end
                try
                    meta.binMapping = strcmp(get_att('binMapping_applied'), 'Bin mapping applied');
                    binmapped = logical(meta.binMapping);
                catch
                    binmapped = false;
                end
                
               
                %subset the global metadata fields to only the respective group.
                dataset_meta_id = ['_' meta_attr_midname '_'];
                [~, other_datasets_meta_names] = filterFields(file_metadata, dataset_meta_id);
                dataset_meta = rmfield(file_metadata, other_datasets_meta_names);
                
                %load extra metadata and unify the variable access pattern into
                % the same function name.
                if is_netcdf
                    meta.dim_meta = data_metadata.(group_name).Dimensions;
                    meta.var_meta = data_metadata.(group_name).Variables;
                    gid = dataset_groups(k);
                    get_var = @(x)(nc_get_var(gid, var_mapping.(x)));
                else
                    fname = getindex(dataset_groups, k);
                    get_var = @(x)(transpose(matdata.(fname).(var_mapping.(x))));
                end
                
                meta.featureType = '';
                meta.instrument_make = 'Nortek';
                meta.instrument_model = get_att('instrument_model');
                
                if is_netcdf
                    inst_serial_numbers = get_att('instrument_serial_no');
                    if numel(unique(inst_serial_numbers)) > 1
                        dispmsg('Multi instrument serial numbers found in %s. Assuming the most frequent is the right one...', filename)
                        inst_serial_no = mode(inst_serial_numbers);
                    else
                        inst_serial_no = inst_serial_numbers(1);
                    end
                else
                    %serial no is at metadata/Config level in the mat files.
                    inst_serial_numbers = get_att('instrument_serial_no');
                    if numel(unique(inst_serial_numbers)) > 1
                        dispmsg('Multi instrument serial numbers found in %s. Assuming the most frequent is the right one...', filename)
                        inst_serial_no = mode(inst_serial_numbers);
                    else
                        inst_serial_no = inst_serial_numbers(1);
                    end
                end
                
                meta.instrument_serial_no = num2str(inst_serial_no);
                
                try
                    assert(contains(meta.instrument_model, 'Signature'))
                    %TODO: support other models. need more files.
                catch
                    errormsg('Only Signature ADCPs are supported.', filename)
                end
                
                default_beam_angle = OceanContourWaves.beam_angles.(meta.instrument_model);
                instrument_beam_angles = single(get_att('beam_angle'));
                try
                    %the attribute may contain 5 beams (e.g. wave).
                    dataset_beam_angles = instrument_beam_angles(1:meta.nBeams);
                    assert(isequal(unique(dataset_beam_angles), default_beam_angle))
                    %TODO: workaround for inconsistent beam_angles. need more files.
                catch
                    errormsg('Inconsistent beam angle/Instrument information in %s', filename)
                end
                meta.beam_angle = default_beam_angle;
                
                meta.('instrument_sample_interval') = single(get_att('instrument_sample_interval'));
                
                %mode_sampling_duration_str = ['instrument_' meta_attr_midname '_interval'];
                %meta.(mode_sampling_duration_str) = get_att(mode_sampling_duration_str);
                % for waves
                meta.('instrument_sampling_duration') = get_att('nSamples') / get_att('sampleRate');
                time = get_var('time')/86400.0 + datenum(1970, 1, 1, 0, 0, 0); %"seconds since 1970-01-01T00:00:00 UTC";
                
                try
                    actual_sample_interval = single(mode(diff(time)) * 86400.);
                    assert(isequal(meta.('instrument_sample_interval'), actual_sample_interval))
                catch
                    expected = meta.('instrument_sample_interval');
                    dispmsg('Inconsistent instrument sampling interval in %s . Metadata is set to %d, while time variable indicates %d. Using variable estimates...', filename, expected, actual_sample_interval);
                    meta.('instrument_sample_interval') = actual_sample_interval;
                end
                
                coordinate_system = get_att('coordinate_system');
                switch coordinate_system
                    case 'XYZ'
                        if logical(get_att('converted_to_enu'))
                            meta.coordinate_system = 'ENU';
                        else
                            errormsg('Unsuported coordinates. %s contains non-ENU data.', filename)
                        end
                    case 'ENU'
                        meta.coordinate_system = 'ENU';
                        % OK
                    otherwise
                        errormsg('Unsuported coordinates. %s contains non-ENU data.', filename)
                end
                
                %                 z = get_var('HEIGHT_ABOVE_SENSOR');
                %                 try
                %                     assert(all(z > 0));
                %                 catch
                %                     errormsg('invalid VelocityENU_Range in %s', filename)
                %                     %TODO: plan any workaround for diff ranges. files!?
                %                 end
                
                %                 binSize = get_var('binSize');
                %                 if numel(unique(binSize)) > 1
                %                     dispmsg('Inconsistent binSizes in %s. Assuming the most frequent is the right one...',filename)
                %                     meta.binSize = mode(binSize);
                %                 else
                %                     meta.binSize = binSize;
                %                 end
                
                meta.file_meta = file_metadata;
                meta.dataset_meta = dataset_meta;
                
                %                 switch meta.coordinate_system
                %                     case 'ENU'
                %                         dimensions = IMOS.gen_dimensions('adcp_enu');
                %                     otherwise
                %                         dimensions = IMOS.gen_dimensions('adcp');
                %                 end
                
                meta.adcp_orientation = 'ZUP';
                
                dataset = struct();
                dataset.toolbox_input_file = filename;
                dataset.toolbox_parser = mfilename;
                dataset.netcdf_group_name = group_name;
                dataset.meta = meta;
                
                % add dimensions with their data
                %wave_dim_names = [wave_dim_names 'LATITUDE' 'LONGITUDE' ];
                nDims = numel(wave_dim_names);
                dataset.dimensions = cell(nDims, 1);
                for i=1:nDims
                    dname = wave_dim_names{i};
                    try
                        imos_name = OceanContourWaves.get_imos_mapped_name(dname, var_map, custom_magnetic_declination);
                    catch
                        imos_name = dname;
                    end
                    dataset.dimensions{i}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(imos_name, 'type')));
                    dataset.dimensions{i}.name         = imos_name;
                    %dataset.dimensions{i}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(imos_name, 'type')));
                    if strcmpi(dname, 'TIME')
                        % TODO: find code that handles cftime conventions
                        % until then assume it will always be
                        % units = "seconds since 1970-01-01T00:00:00 UTC"
                        dataset.dimensions{i}.data = dataset.dimensions{i}.typeCastFunc(get_var(dname))/86400 + datenum(1970,1,1,0,0,0);
                    else
                        dataset.dimensions{i}.data = dataset.dimensions{i}.typeCastFunc(get_var(dname));
                    end
                    %                   if strcmpi(dims{i, 1}, 'DIR')
                    %                       dataset.dimensions{i}.compass_correction_applied = meta.compass_correction_applied;
                    %                       dataset.dimensions{i}.comment  = magdec_attrs.comment;
                    %                   end
                end
                
                % list of variable names (not including known dimension
                % names

                idx = contains(known_var_names, [wave_dim_names {'time'} {'TIME'}]);
                var_names = known_var_names(~idx);
                
                % add variables with their data mapped
                nVars = numel(var_names) + 3;
                dataset.variables = cell(nVars, 1);
                dataset.variables{1}.name = 'TIMESERIES';
                dataset.variables{1}.dimensions = [];
                dataset.variables{1}.data = 1;
                dataset.variables{1}.typeCastFunc =  str2func('int32');
                dataset.variables{2}.name = 'LATITUDE';
                dataset.variables{2}.dimensions = [];
                dataset.variables{2}.data = NaN;
                dataset.variables{2}.typeCastFunc =  str2func('double');
                dataset.variables{3}.name = 'LONGITUDE';
                dataset.variables{3}.dimensions = [];
                dataset.variables{3}.data = NaN;
                dataset.variables{3}.typeCastFunc =  str2func('double');
                dataset.variables{4}.name = 'NOMINAL_DEPTH';
                dataset.variables{4}.dimensions = [];
                dataset.variables{4}.data = NaN;
                dataset.variables{4}.typeCastFunc =  str2func('double');
                
                for i=4:nVars

                    vname = var_names{i-3};
                    idx = strcmp(vname, known_var_names);
                    vstruct = json_varwaves.variables(idx);
                    
                    if isstruct(vstruct.attribute)
                        vstruct.attribute = num2cell(vstruct.attribute);
                    end
                    imos_name = OceanContourWaves.get_imos_mapped_name(vname, var_map, custom_magnetic_declination);
                    dataset.variables{i}.name = imos_name;
                    dataset.variables{i}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(imos_name, 'type')));
                    % get dimension names of the variable and map to IMOS
                    % names
                    vdimnames = strsplit(vstruct.shape, ' ');
                    vdimnames_imos = vdimnames;
                    for j = 1:numel(vdimnames)
                        dname = vdimnames{j};
                        imos_dname = OceanContourWaves.get_imos_mapped_name(dname, var_map, custom_magnetic_declination);
                        vdimnames_imos{j} = imos_dname;
                    end
                    % get indicies into dimension array
                    dim_ind = cellfun(@(x) find(ismember(wave_dim_names, x)), vdimnames);
                    dataset.variables{i}.dimensions = dim_ind;
                    % some variables have dims [adimname time]
                    
                    if numel(dim_ind) == 1
                        dataset.variables{i}.data = dataset.variables{i}.typeCastFunc(get_var(vname));
                    else
                        % Some variables have dimension order list with time at the end eg "EnergySpectra_Frequency time"
                        % I don't know if this also would happen in mat
                        % export
                        if dim_ind(end) == 1
                            dataset.variables{i}.dimensions = fliplr(dim_ind);
                            dataset.variables{i}.data = dataset.variables{i}.typeCastFunc(get_var(vname));
                        else
                            dataset.variables{i}.data = permute(dataset.variables{i}.typeCastFunc(get_var(vname)), numel(dim_ind):-1:1);
                        end
                    end
                    dataset.variables{i}.coordinates = 'TIME LATITUDE LONGITUDE';
                    % get variable description and save as comment
                    attr_names = cellfun(@(x) x.name, vstruct.attribute, 'UniformOutput', false);
                    [tf, ind] = inCell(attr_names, 'description');
                    if tf
                        attr = vstruct.attribute{ind};
                        if ~isempty(attr.value)
                            dataset.variables{i}.comment = attr.value;
                        end
                    end
                    % get units and save if valid
                    [tf, ind] = inCell(attr_names, 'units');
                    if tf
                        attr = vstruct.attribute{ind};
                        if ~isempty(attr.value) && ~strcmp(attr.value, '?')
                            dataset.variables{i}.units = attr.value;
                        end
                    end
                    %                   if strcmpi(dims{i, 1}, 'DIR')
                    %                       dataset.dimensions{i}.compass_correction_applied = meta.compass_correction_applied;
                    %                       dataset.dimensions{i}.comment  = magdec_attrs.comment;
                    %                   end
                end
                
                sample_data{k} = dataset;
            end
            
            if is_netcdf
                netcdf.close(ncid);
            end
        end
        
    end
    
end
