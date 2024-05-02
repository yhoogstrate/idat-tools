#!/usr/bin/env python

# get sections
# https://github.com/FoxoTech/methylprep/blob/master/methylprep/files/idat.py#L228

# 1. open file

#get_file_object() -> idat_file
#open(filepath_or_buffer, 'rb')


from parser import read_char, read_long


def get_magic(fh_in):
    fh_in.seek(0) # IdatHeaderLocation.FILE_TYPE.value
    file_type = read_char(fh_in, 4)
    
    return file_type

def get_idat_version(fh_in):
    fh_in.seek(4) # IdatHeaderLocation.VERSION.value
    idat_version = read_long(fh_in)
    
    return idat_version


with open("207513420127_R08C01_Grn.idat", "rb") as fh_in:
    print("1.  Magic:        ["+get_magic(fh_in)+"]")
    print("2.  IDAT version: ["+str(get_idat_version(fh_in))+"]")
    

