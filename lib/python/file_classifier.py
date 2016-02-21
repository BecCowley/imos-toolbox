"""FileClassifier - Generic class for working out the destination path of
a file to be published. The idea is to define the common functionality
here, then create subclasses to customise for each specific incoming
handler. 

Expected use:

    class MyFileClassifier(FileClassifier):
        def dest_path(self, input_file):
            path = <case-specific logic> 
            ...
            return path

    try:
        dest_path = MyFileClassifier.dest_path(input_file)
    except FileClassifierException, e:
        print >>sys.stderr, e
        exit(1)

    print dest_path

"""

import os
import sys
import re
from netCDF4 import Dataset


class FileClassifierException(Exception):
    pass

class FileClassifier(object):
    "Base class for working out where a file should be published."

    @classmethod
    def _error(cls, message):
        "Raise an exception with the given message."
        raise FileClassifierException, message

    @classmethod
    def _open_nc_file(cls, file_path):
        "Open a NetCDF file for reading"
        try:
            return Dataset(file_path, mode='r')
        except:
            cls._error("Could not open NetCDF file '%s'." % file_path)

    @classmethod
    def _get_nc_att(cls, file_path, att_name, default=None):
        """Return the value of a global attribute from a NetCDF file. If a
        list of attribute names is given, a list of values is
        returned.  Unless a default value other than None is given, a
        missing attribute raises an exception.

        """
        dataset = cls._open_nc_file(file_path)

        if isinstance(att_name, list):
            att_list = att_name
        else:
            att_list = [att_name]
        values = []

        for att in att_list:
            val = getattr(dataset, att, default)
            if val is None:
                cls._error("File '%s' has no attribute '%s'" % (file_path, att))
            values.append(val)
        dataset.close()

        if isinstance(att_name, list):
            return values
        return values[0]


    @classmethod
    def _get_site_code(cls, input_file):
        "Return the site_code attribute of the input_file"
        return cls._get_nc_att(input_file, 'site_code')

    @classmethod
    def _make_path(cls, dir_list):
        """Create a path from a list of directory names, making sure the
         result is a plain ascii string, not unicode (which could
         happen if some of the components of dir_list come from NetCDF
         file attributes).

        """
        for i in range(len(dir_list)):
            dir_list[i] = str(dir_list[i])
        return os.path.join(*dir_list)