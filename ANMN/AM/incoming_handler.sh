#!/bin/bash

export PYTHONPATH="$DATA_SERVICES_DIR/ANMN"
export SCRIPTPATH="$DATA_SERVICES_DIR/ANMN/AM"

# returns extension of file
# $1 - file
_get_extension() {
    local file=$1; shift
    echo ${file##*.}
}

# handle a netcdf file for the facility
# $1 - file to handle
handle_netcdf() {
    local file=$1; shift

    check_netcdf          $file         || file_error $file "Not a valid NetCDF file"
    check_netcdf_cf       $file         || file_error $file "NetCDF file is not CF compliant"
    check_netcdf_imos     $file         || file_error $file "NetCDF file is not IMOS compliant"
    check_netcdf_facility $file anmn_am || file_error $file "NetCDF file is not ANMN_AM compliant"

    local path_hierarchy
    path_hierarchy=`$SCRIPTPATH/destPath.py $file` || file_error $file "Could not determine destination path for file"
    [ x"$path_hierarchy" = x ] && file_error $file "Could not determine destination path for file"

    # TODO: archive previous version of file if found on opendap

    move_to_opendap_imos $file $path_hierarchy
}

# handle a netcdf file for the facility
# $1 - file to handle
handle_csv() {
    local file=$1; shift

    # generate NetCDF file using python script
    local netcdf_file
    local wip_dir="$WIP_DIR/ANMN/AM"
    mkdir -p $wip_dir || file_error "Could not create wip directory '$wip_dir'"
    netcdf_file=`cd $wip_dir && $SCRIPTPATH/rtCO2.py $file` || file_error $file "Could not generate NetCDF file"

    if [ x"$netcdf_file" = x ]; then
        # no new NetCDF file created because csv file contains no new data since last run
        log_info "Nothing new to process"
        rm -f $file
    else
        local netcdf_file_full_path="$wip_dir/$netcdf_file"
        test -f $netcdf_file_full_path || file_error $file "Could not generate NetCDF file"

        handle_netcdf $netcdf_file_full_path && \
            rm -f $file
    fi
}

# main
# $1 - file to handle
main() {
    local file=$1; shift

    local file_extension=`_get_extension $file`

    if [ "$file_extension" = "nc" ]; then
        handle_netcdf $file
    elif [ "$file_extension" = "csv" ]; then
        handle_csv $file
    else
        file_error $file "Not a NetCDF nor a csv file, extension is '$file_extension'"
    fi
}


main "$@"
