classdef OceanContour
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
            % midname = OceanContour.build_meta_attr_midname('Avg');
            % assert(strcmp(midname,'avg'))
            % midname = OceanContour.build_meta_attr_midname('burstAltimeter');
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
            % name = OceanContour.build_instrument_name('Avg', 'coordSystem');
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

            meta_attr_midname = OceanContour.build_meta_attr_midname(group_name);
            attname = ['Instrument_' meta_attr_midname '_' var_name];
        end

        function [ucur_name, vcur_name, heading_name] = build_magnetic_variables(custom_magnetic_declination)
            %function attname = build_magnetic_variables(custom_magnetic_declination)
            %
            % Generate VAR or VAR_MAG toolbox variable style names
            % based on provided magnetic declination info.
            %
            narginchk(1, 1)

            if ~islogical(custom_magnetic_declination)
                errormsg('build_magnetic_variables: first argument is not a logical')
            end

            if ~custom_magnetic_declination
                %TODO: This is probably unecessary
                %I believe OceanContourDouble-check if OceanContour will change variable names if custom magnetic declination is used.
                dispmsg('%s: Assigning non ENU Velocities to ENU variables. Verify the magnetic declination angles.')
                ucur_name = 'UCUR_MAG';
                vcur_name = 'VCUR_MAG';
                heading_name = 'HEADING_MAG';
            else
                ucur_name = 'UCUR';
                vcur_name = 'VCUR';
                heading_name = 'HEADING';
            end

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
            % attmap = OceanContour.get_attmap('Avg');
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

            meta_attr_midname = OceanContour.build_meta_attr_midname(group_name);

            % attribute names seem to be different for different files, eg,
            % CSIRO ones don't have the 'Instrument' prefix for all. Let's
            % try finding the suffixes and matching that way

            flds = fields(file_metadata);
            iignore = endsWith(flds,'description');
            descflds = flds(iignore);
            flds = flds(~iignore);
            attmap.('instrument_model') = flds{contains(flds,'instrumentName')};
            attmap.('beam_angle') = flds{contains(flds,'slantAngles')};
            %attmap.('beam_interval') = OceanContour.build_instrument_name(group_name, 'slantAngles');
            attmap.('coordinate_system') = flds{contains(flds,'coordSystem')};
            attmap.('converted_to_enu') = descflds{contains(descflds,'transformsAndCorrections_addENU_description')};
            
            % while this might be a waves file, some info is 'avg' prefix
            if isfield(file_metadata, 'Instrument_avg_enable') & logical(file_metadata.Instrument_avg_enable)
                inst_data_name = 'avg';
            else
                inst_data_name = 'burst';
            end
            has_data_been_averaged = false;
            iaverage = find(contains(flds,'average_data'));
            if ~isempty(iaverage) 
                if logical(file_metadata.(flds{iaverage}))
                    has_data_been_averaged = true;
                end
            end
            
            attmap.('nBeams') = flds{contains(flds,'nBeams')};
            attmap.('activeBeams') = flds{contains(flds,'activeBeams')}; %no previous name
            attmap.('magDec_DataInfo') = flds{contains(flds,'trig_en')};
            attmap.('magDec_User') = flds{contains(flds,'Instrument_user_decl')};
            attmap.('binMapping') = flds{contains(flds,'transformsAndCorrections_binMapping')};            
            attmap.('binMapping_applied') = descflds{contains(descflds,'transformsAndCorrections_binMapping_description')};    
 
            if strcmpi(ftype, 'mat')
                attmap.('instrument_serial_no') = 'Instrument_serialNumberDoppler';
                attmap.('binSize') = OceanContour.build_instrument_name(group_name, 'cellSize');
            else
                attmap.('instrument_serial_no') = flds{contains(flds,'serial')};
                attmap.('binSize') = flds{contains(flds,['Instrument_' meta_attr_midname '_cellSize'])};

            end

            %custom & dynamical fields
            attmap.(['instrument_' meta_attr_midname '_enable']) = OceanContour.build_instrument_name(group_name, 'enable');

            switch meta_attr_midname
                case 'avg'
                    attmap.('instrument_avg_interval') = OceanContour.build_instrument_name(group_name, 'averagingInterval');
                    attmap.('instrument_sample_interval') = OceanContour.build_instrument_name(group_name, 'measurementInterval');
                    %TODO: need a more complete file to test below below
                case 'burst'
                    % TODO: if burst data that has been averaged, what is
                    % the timestep
                    % If DataInfo_average_data set and is delta time = DataInfo_average_window
                    % or delta time = DataInfo_average_window/2?
                    attmap.('instrument_burst_interval') = OceanContour.build_instrument_name(group_name, 'measurementInterval');
                    attmap.('instrument_sample_rate') = flds{contains(flds, 'sampleRate')};
                    if has_data_been_averaged
                        attmap.('instrument_sample_interval') = flds{contains(flds,'average_window')};
                    end
                case 'bursthr'
                    attmap.('instrument_bursthr_interval') = OceanContour.build_instrument_name(group_name, 'burstHourlyInterval');
                case 'burstAltimeter'
                    attmap.('instrument_burstAltimeter_interval') = OceanContour.build_instrument_name(group_name, 'burstAltimeterInterval');
                case 'burstRawAltimeter'
                    attmap.('instrument_burstRawAltimeter_interval') = OceanContour.build_instrument_name(group_name, 'burstRawAltimeterInterval');
            end

        end

        function [varmap] = get_varmap(ftype, vars, nbeams, custom_magnetic_declination, binmapped, is_enu)
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
            % varmap = OceanContour.get_attmap('Avg',4,False);
            % assert(strcmp(attmap.WCUR_2,'Vel_Up2'));
            %
            % % nbeams == 3
            % varmap = OceanContour.get_varmap('Avg',3,False);
            % f=false;try;varmap.WCUR_2;catch;f=true;end
            % assert(f)
            %
            % % custom magdec - may change with further testing
            % varmap = OceanContour.get_varmap('Avg',4,True);
            % assert(strcmp(varmap.UCUR_MAG,'Vel_East'))
            %
            %
            % author: hugo.oliveira@utas.edu.au
            %
            narginchk(6, 6)

            if ~ischar(ftype)
                errormsg('First argument is not a string')
            elseif ~strcmpi(ftype, 'mat') && ~strcmpi(ftype, 'netcdf')
                errormsg('First argument %s is an invalid ftype. Accepted file types are ''mat'' and ''netcdf''.', ftype)
            elseif ~isstruct(vars)
                errormsg('Second argument is not a structure')
            elseif ~isscalar(nbeams)
                errormsg('Third argument is not a scalar')
            elseif ~islogical(custom_magnetic_declination)
                errormsg('Fourth argument is not logical')
            elseif ~islogical(binmapped)
                errormsg('Fifth argument is not logical')
            end

            is_netcdf = strcmpi(ftype, 'netcdf');
            [ucur_name, vcur_name, heading_name] = OceanContour.build_magnetic_variables(custom_magnetic_declination);
            
            flds = fieldnames(vars);
            varmap = struct();
            varmap.('binSize') = 'CellSize';
            varmap.('TIME') = 'MatlabTimeStamp';
            varmap.('TIME_CFTIME') = 'time';
            
            if is_netcdf
                varmap.('status') = 'Status';
                %TODO: reinforce uppercase at first letter? nEed to see more files.
                varmap.('HEIGHT_ABOVE_SENSOR') = flds{contains(flds, 'VelocityENU_Range')};
                %TODO: Handle magnetic & along beam cases.
                %varmap.('DIST_ALONG_BEAMS') = [group_name 'Velocity???_Range'];
                %TODO: evaluate if when magnetic declination is provided, the
                %velocity fields will be corrected or not (as well as any rename/comments added).
                varmap.(heading_name) = 'Heading';
                if is_enu
                    varmap.(ucur_name) = flds{contains(flds, 'Vel_East')};
                    varmap.(vcur_name) = flds{contains(flds, 'Vel_North')};
                    varmap.('WCUR') = flds{contains(flds, 'Vel_Up1')};
                else
                    varmap.('VEL1') = 'Vel_Beam1';
                    varmap.('VEL2') = 'Vel_Beam2';
                    varmap.('VEL3') = 'Vel_Beam3';                    
                end
                varmap.('ABSI1') = flds{contains(flds, 'Amp_Beam1')};
                varmap.('ABSI2') = flds{contains(flds,'Amp_Beam2')};
                varmap.('ABSI3') = flds{contains(flds, 'Amp_Beam3')};
                varmap.('CMAG1') = flds{contains(flds, 'Cor_Beam1')};
                varmap.('CMAG2') = flds{contains(flds, 'Cor_Beam2')};
                varmap.('CMAG3') = flds{contains(flds, 'Cor_Beam3')};

                if nbeams > 3
                    if is_enu
                        varmap.('WCUR_2') = flds{contains(flds,'Vel_Up2')};
                    else
                        varmap.('VEL4') = 'Vel_Beam4';
                    end
                    varmap.('ABSI4') = flds{contains(flds, 'Amp_Beam4')};
                    varmap.('CMAG4') = flds{contains(flds, 'Cor_Beam4')};
                end

            else
                %TODO: check if norteks also change the variable names
                %when exporting to matlab.
                %instrument_serial_no is on metadata for matfiles.
                varmap.('status') = 'Status';
                varmap.('HEIGHT_ABOVE_SENSOR') = 'Range';
                varmap.(heading_name) = 'Heading';
                if is_enu
                    varmap.(ucur_name) = flds{contains(flds, 'VelEast')};
                    varmap.(vcur_name) = flds{contains(flds, 'VelNorth')};
                    varmap.('WCUR') = flds{contains(flds, 'VelUp1')};
                else
                    varmap.('VEL1') = 'VelBeam1';
                    varmap.('VEL2') = 'VelBeam2';
                    varmap.('VEL3') = 'VelBeam3';
                end
                varmap.('ABSI1') = flds{contains(flds, 'AmpBeam1')};
                varmap.('ABSI2') = flds{contains(flds, 'AmpBeam2')};
                varmap.('ABSI3') = flds{contains(flds, 'AmpBeam3')};
                varmap.('CMAG1') = flds{contains(flds, 'CorBeam1')};
                varmap.('CMAG2') = flds{contains(flds, 'CorBeam2')};
                varmap.('CMAG3') = flds{contains(flds, 'CorBeam3')};

                if nbeams > 3
                    if is_enu
                        varmap.('WCUR_2') = flds{contains(flds, 'VelUp2')};
                    else
                        varmap.('VEL4') = 'VelBeam4';
                    end
                    varmap.('ABSI4') = flds{contains(flds, 'AmpBeam4')};
                    varmap.('CMAG4') = flds{contains(flds, 'CorBeam4')};
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

        function [imap] = get_importmap(nbeams, custom_magnetic_declination)
            %function [imap] = get_importmap(custom_magnetic_declination)
            %
            % Return default variables to import from the OceanContour files.
            %
            % Inputs:
            %
            % nbeams [scalar] - the number of ADCP beams.
            % custom_magnetic_declination [logical] - true for custom
            %                                         magnetic values.
            %
            % Outputs:
            %
            % imap [struct[cell]] - Struct with different variables
            %                       classes to import
            %
            %
            % Example:
            %
            % %basic usage
            % imap = OceanContour.get_importmap(False);
            % assert(inCell(imap.all_variables,'PITCH'))
            % assert(inCell(imap.all_variables,'ROLL'))
            %
            % author: hugo.oliveira@utas.edu.au
            %
            narginchk(2, 2)

            if ~isscalar(nbeams)
                errormsg('First argument is not a scalar')
            elseif ~islogical(custom_magnetic_declination)
                errormsg('Second argument is not a logical')
            end

            imap = struct();
            [ucur_name, vcur_name, heading_name] = OceanContour.build_magnetic_variables(custom_magnetic_declination);

            ENU = struct();

            ENU.one_dimensional = {'TEMP', 'PRES_REL', 'SSPD', 'BAT_VOLT', 'PITCH', 'ROLL', heading_name, 'ERROR', 'AMBIG_VEL', 'TRANSMIT_E', 'NOMINAL_CORR'};
            ENU.velocity_variables = {ucur_name, vcur_name, 'WCUR'};
            ENU.beam_amplitude_variables = {'ABSI1', 'ABSI2', 'ABSI3'};
            ENU.correlation_variables = {'CMAG1', 'CMAG2', 'CMAG3'};

            if nbeams > 3
                ENU.velocity_variables = [ENU.velocity_variables, 'WCUR_2'];
                ENU.beam_amplitude_variables = [ENU.beam_amplitude_variables 'ABSI4'];
                ENU.correlation_variables = [ENU.correlation_variables 'CMAG4'];
            end

            ENU.two_dimensional = [ENU.velocity_variables, ENU.beam_amplitude_variables, ENU.correlation_variables];
            ENU.all_variables = [ENU.one_dimensional, ENU.two_dimensional];

            %TODO: Implement Non-ENU cases.

            imap.('ENU') = ENU;

            BEAM = struct();

            BEAM.one_dimensional = {'TEMP', 'PRES_REL', 'SSPD', 'BAT_VOLT', 'PITCH', 'ROLL', heading_name, 'ERROR', 'AMBIG_VEL', 'TRANSMIT_E', 'NOMINAL_CORR'};
            BEAM.velocity_variables = {'VEL1', 'VEL2', 'VEL3'};
            BEAM.beam_amplitude_variables = {'ABSI1', 'ABSI2', 'ABSI3'};
            BEAM.correlation_variables = {'CMAG1', 'CMAG2', 'CMAG3'};

            if nbeams > 3
                BEAM.velocity_variables = [BEAM.velocity_variables, 'VEL4'];
                BEAM.beam_amplitude_variables = [BEAM.beam_amplitude_variables 'ABSI4'];
                BEAM.correlation_variables = [BEAM.correlation_variables 'CMAG4'];
            end

            BEAM.two_dimensional = [BEAM.velocity_variables, BEAM.beam_amplitude_variables, BEAM.correlation_variables];
            BEAM.all_variables = [BEAM.one_dimensional, BEAM.two_dimensional];
            
            imap.('BEAM') = BEAM;
            
            XYZ = struct();

            XYZ.one_dimensional = {'TEMP', 'PRES_REL', 'SSPD', 'BAT_VOLT', 'PITCH', 'ROLL', heading_name, 'ERROR', 'AMBIG_VEL', 'TRANSMIT_E', 'NOMINAL_CORR'};
            XYZ.velocity_variables = {'VEL1', 'VEL2', 'VEL3'};
            XYZ.beam_amplitude_variables = {'ABSI1', 'ABSI2', 'ABSI3'};
            XYZ.correlation_variables = {'CMAG1', 'CMAG2', 'CMAG3'};

            if nbeams > 3
                XYZ.velocity_variables = [XYZ.velocity_variables, 'VEL4'];
                XYZ.beam_amplitude_variables = [XYZ.beam_amplitude_variables 'ABSI4'];
                XYZ.correlation_variables = [XYZ.correlation_variables 'CMAG4'];
            end

            XYZ.two_dimensional = [XYZ.velocity_variables, XYZ.beam_amplitude_variables, XYZ.correlation_variables];
            XYZ.all_variables = [XYZ.one_dimensional, XYZ.two_dimensional];
            
            imap.('XYZ') = XYZ;
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
            % [sample_data] = readOceanContour(file);
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
            % [sample_data] = readOceanContour(file);
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
                try
                    matdata = load(filename);
                    ftype = 'mat';
                catch
                    errormsg('%s is not a mat or netcdf file', filename)
                end

            end

            is_netcdf = strcmpi(ftype, 'netcdf');

            if is_netcdf
                OceanContour.verify_netcdf_groups(info);
                file_metadata = nc_flat(info.Groups(1).Attributes, false);
                data_metadata = nc_flat(info.Groups(2).Groups, false);

                ncid = netcdf.open(filename);
                c = onCleanup(@()netcdf.close(ncid));
                root_groups = netcdf.inqGrps(ncid);
                data_group = root_groups(2);

                dataset_groups = netcdf.inqGrps(data_group);
                get_group_name = @(x)(netcdf.inqGrpName(x));
                
            else
                OceanContour.verify_mat_groups(matdata);
                file_metadata = matdata.Config;
                matdata = rmfield(matdata, 'Config'); %mem optimisation.

                dataset_groups = fieldnames(matdata);
                get_group_name = @(x)(getindex(split(x, '_Data'), 1));

            end

            n_datasets = numel(dataset_groups);
            sample_data = cell(1, n_datasets);

            for k = 1:n_datasets

                % start by loading preliminary information into the metadata struct, so we
                % can define the variable names and variables to import.
                meta = struct();

                group_name = get_group_name(dataset_groups);
                meta_attr_midname = OceanContour.build_meta_attr_midname(group_name);

                %load toolbox_attr_names:file_attr_names dict.
                att_mapping = OceanContour.get_attmap(file_metadata, ftype, group_name);

                %access pattern - use lookup based on expected names,
                get_att = @(x)(file_metadata.(att_mapping.(x)));

                nBeams = double(get_att('nBeams'));

                try
                    activeBeams = double(get_att('activeBeams'));
                catch
                    activeBeams = Inf;
                end

                meta.nBeams = min(nBeams, activeBeams);

                try
                    assert(meta.nBeams == 4 | meta.nBeams == 5);
                    %TODO: support variable nBeams. need more files.
                catch
                    errormsg('Only 4 Beam ADCP are supported. %s got %d nBeams', filename, meta.nBeams)
                end

                magDec_User = get_att('magDec_User');
                has_magdec_user = logical(magDec_User);
                try
                    magDec_DataInfo = get_att('magDec_DataInfo');
                    has_magdec_oceancontour = logical(magDec_DataInfo);
                catch
                    magDec_DataInfo = NaN;
                    has_magdec_oceancontour = false;
                end
                meta.magDec = 0.0;
                custom_magnetic_declination = has_magdec_user | has_magdec_oceancontour;
                if has_magdec_oceancontour
                    meta.magDec = magDec_DataInfo;
                elseif has_magdec_user
                    meta.magDec = magDec_User;
                end
                
                try
                    meta.binMapping = get_att('binMapping');
                catch
                    meta.binMapping = false;
                end

                try
                    meta.binMapping = strcmp(get_att('binMapping_applied'), 'Bin mapping applied');
                    binmapped = logical(meta.binMapping);
                catch
                    binmapped = false;
                end
                
                is_waves = false;
                if isfield(file_metadata, 'DataInfo_waves_processing')
                    is_waves = true;
                end
                
                %Now that we know some preliminary info, we can load the variable
                % name mappings and the list of variables to import.

                coordinate_system = get_att('coordinate_system');
                switch coordinate_system
                    case {'XYZ', 'BEAM'}
                        meta.coordinate_system = coordinate_system;
                        try
                            has_converted_to_enu = logical(get_att('converted_to_enu'));
                        catch
                            has_converted_to_enu = false;
                        end
                        if has_converted_to_enu
                            meta.coordinate_system = 'ENU';
                        else
                            dispmsg('Unsuported coordinates. %s contains non-ENU data.', filename)
                        end
                    case 'ENU'
                        meta.coordinate_system = 'ENU';
                        % OK
                    otherwise
                        errormsg('Unsuported coordinates. %s contains non-ENU data.', filename)
                end
                is_enu = strcmp(meta.coordinate_system, 'ENU');
                var_mapping = OceanContour.get_varmap(ftype, group_name, nBeams, custom_magnetic_declination, is_binmapped, is_enu);
                import_mapping = OceanContour.get_importmap(nBeams, custom_magnetic_declination);

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

                % Have observed occasional unfinished/invalid
                % packet/ensemble. These invalid packets so far always have
                % a status == 0. So store 
                iGood = logical(get_var('status'));
				
                if is_netcdf
                    inst_serial_numbers = get_var('instrument_serial_no');
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

                default_beam_angle = OceanContour.beam_angles.(meta.instrument_model);
                instrument_beam_angles = single(get_att('beam_angle'));
                try
                    %the attribute may contain 5 beams (e.g. AST for waves).
                    % TODO: workaround for inconsistent beam_angles. need more files.
                    % At the moment just assume the first nBeams-1 are for 
                    % velocity calculation
                    dataset_beam_angles = instrument_beam_angles(1:meta.nBeams);
                    assert(isequal(unique(dataset_beam_angles(1:meta.nBeams-1)), default_beam_angle))
                    %TODO: workaround for inconsistent beam_angles. need more files.
                catch
                    errormsg('Inconsistent beam angle/Instrument information in %s', filename)
                end
                meta.beam_angle = default_beam_angle;

                meta.('instrument_sample_interval') = single(get_att('instrument_sample_interval')/get_att('instrument_sample_rate'));

                % TODO
                %mode_sampling_duration_str = ['instrument_' meta_attr_midname '_interval'];
                %meta.(mode_sampling_duration_str) = get_att(mode_sampling_duration_str);

                time = get_var('TIME');
                time_cftime = nc_get_var(gid, 'time')/86400.0 + datenum(1970,1,1,0,0,0); %"seconds since 1970-01-01T00:00:00 UTC";
                
                try
                    actual_sample_interval = single(mode(diff(time)) * 86400.);
                    assert(isequal(meta.('instrument_sample_interval'), actual_sample_interval))
                catch
                    expected = meta.('instrument_sample_interval');                    
                    dispmsg('Inconsistent instrument sampling interval in %s . Metadata is set to %d, while time variable indicates %d. Using variable estimates...', filename, expected, actual_sample_interval);
                    meta.('instrument_sample_interval') = actual_sample_interval;                    
                end
              
                z = get_var('HEIGHT_ABOVE_SENSOR');              
                try
                    assert(all(z > 0));
                catch
                    errormsg('invalid VelocityENU_Range in %s', filename)
                    %TODO: plan any workaround for diff ranges. files!?
                end

                binSize = get_var('binSize');                
                if numel(unique(binSize)) > 1
                    dispmsg('Inconsistent binSizes in %s. Assuming the most frequent is the right one...',filename)                    
                    meta.binSize = mode(binSize);                    
                else                    
                    meta.binSize = binSize;
                end                                   

                meta.file_meta = file_metadata;
                meta.dataset_meta = dataset_meta;
                                
                switch meta.coordinate_system
                    case 'ENU'                
                        dimensions = IMOS.gen_dimensions('adcp_enu');                                                        
                    otherwise                        
                        dimensions = IMOS.gen_dimensions('adcp');
                end
                
                status_data = get_var('status');
                adcpOrientations = arrayfun(@(x) bin2dec(num2str(bitget(x, 28:-1:26, 'uint32'))), status_data);
                adcpOrientation = mode(adcpOrientations); % hopefully the most frequent value reflects the orientation when deployed
                % we assume adcpOrientation == 4 by default "ZUP"
                meta.adcp_orientation = 'ZUP';
                adcp_orientation_conversion  = 1;
                if adcpOrientation == 5
                    meta.adcp_orientation = 'ZDOWN';
                    adcp_orientation_conversion  = -1;
                end
                
                dimensions{1}.data = time_cftime;
                dimensions{1}.comment = 'time imported from matlabTimeStamp variable';
                dimensions{2}.data = z * adcp_orientation_conversion ;
                dimensions{2}.comment = 'height imported from VelocityENU_Range';

                switch meta.coordinate_system
                    case 'ENU'
                        onedim_vnames = import_mapping.('ENU').one_dimensional;
                        twodim_vnames = import_mapping.('ENU').two_dimensional;
                    otherwise
                        onedim_vnames = import_mapping.(meta.coordinate_system).one_dimensional;
                        twodim_vnames = import_mapping.(meta.coordinate_system).two_dimensional;
                        dispmsg('%s coordinates found in %s is not implemented yet', filename, meta.coordinate_system)
                end

                onedim_vcoords = [dimensions{1}.name ' LATITUDE LONGITUDE ' 'NOMINAL_DEPTH']; %TODO: point to Pressure/Depth via CF-conventions
                onedim_vtypes = IMOS.cellfun(@getIMOSType, onedim_vnames);
                [onedim_vdata, failed_items] = IMOS.cellfun(get_var, onedim_vnames);

                if ~isempty(failed_items)
                    OceanContour.warning_failed(failed_items, filename)
                end

                twodim_vcoords = [dimensions{1}.name ' LATITUDE LONGITUDE ' dimensions{2}.name];
                twodim_vtypes = IMOS.cellfun(@getIMOSType, twodim_vnames);
                [twodim_vdata, failed_items] = IMOS.cellfun(get_var, twodim_vnames);
                if ~isempty(failed_items)
                    OceanContour.warning_failed(failed_items, filename)
                end
                
                try
                    twodim_vdatamask = get_var('data_mask');
                catch
                    twodim_vdatamask = [];
                end
                has_data_mask = ~isempty(twodim_vdatamask);
                meta.twodim_vdatamask = twodim_vdatamask;
                
                %TODO: Implement unit conversions monads.
                variables = [...
                            IMOS.featuretype_variables('timeSeries'), ...
                            IMOS.gen_variables(dimensions, onedim_vnames, onedim_vtypes, onedim_vdata, 'coordinates', onedim_vcoords), ...
                            IMOS.gen_variables(dimensions, twodim_vnames, twodim_vtypes, twodim_vdata, 'coordinates', twodim_vcoords), ...
                            ];

                dataset = struct();
                dataset.toolbox_input_file = filename;
                dataset.toolbox_parser = mfilename;
                dataset.netcdf_group_name = group_name;
                dataset.meta = meta;
                dataset.dimensions = dimensions;
                dataset.variables = variables;

                sample_data{k} = dataset;
            end

        end

    end

end
