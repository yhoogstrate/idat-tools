#!/usr/bin/env python

# get sections
# https://github.com/FoxoTech/methylprep/blob/master/methylprep/files/idat.py#L228

# 1. open file

#get_file_object() -> idat_file
#open(filepath_or_buffer, 'rb')


def read_char(infile, num_bytes):
    """Converts an array of bytes to a string.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.
        num_bytes {integer} -- The number of bytes to read and parse.

    Returns:
        [string] -- UTF-8 decoded string value.
    """
    return infile.read(num_bytes).decode('utf-8')


def get_magic(fh_in):
    fh_in.seek(0) # IdatHeaderLocation.FILE_TYPE.value
    file_type = read_char(fh_in, 4)
    
    return file_type


"""
def is_idat_file(fh_in):
    idat_file.seek(0) # IdatHeaderLocation.FILE_TYPE.value
    file_type = read_char(idat_file, len(expected))
    return file_type.lower() == expected.lower()
"""


with open("207513420127_R08C01_Grn.idat", "rb") as fh_in:
    print("1. Magic: ["+get_magic(fh_in)+"]")

