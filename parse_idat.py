#!/usr/bin/env python

# get sections
# https://github.com/FoxoTech/methylprep/blob/master/methylprep/files/idat.py#L228

# 1. open file

#get_file_object() -> fh_in
#open(filepath_or_buffer, 'rb')


from parser import *


sections = {
    102: 'ILLUMINA_ID',
    103: 'STD_DEV',
    104: 'MEAN',
    107: 'NUM_BEADS / rep measurements per probe',
    200: 'MID_BLOCK',
    300: 'RUN_INFO',
    400: 'RED_GREEN',
    401: 'MOSTLY_NULL / manifest',
    402: 'BARCODE',
    403: 'CHIP_TYPE / format',
    404: 'MOSTLY_A / label',
    405: 'UNKNOWN_1 / OPA',
    406: 'UNKNOWN_2 / SID',
    407: 'UNKNOWN_3 / DESCR',
    408: 'UNKNOWN_4 / plate ID',
    409: 'UNKNOWN_5 / well ID',
    410: 'UNKNOWN_6',
    510: 'UNKNOWN_7',
    1000: 'NUM_SNPS_READ'
}


def get_magic(fh_in):
    fh_in.seek(0) # IdatHeaderLocation.FILE_TYPE.value
    file_type = read_char(fh_in, 4)
    
    return file_type

def get_idat_version(fh_in):
    fh_in.seek(4) # IdatHeaderLocation.VERSION.value
    idat_version = read_long(fh_in)
    
    return idat_version


def get_section_offsets(fh_in):
    """Parses the IDAT file header to get the byte position
    for the start of each section.

    Arguments:
        fh_in {file-like} -- the IDAT file to process.

    Returns:
        [dict] -- The byte offset for each file section.
    """
    fh_in.seek(12) # IdatHeaderLocation.FIELD_COUNT.value
    num_fields = read_int(fh_in)

    fh_in.seek(16) # IdatHeaderLocation.SECTION_OFFSETS.value

    offsets = {}
    for _idx in range(num_fields):
        key = read_short(fh_in)
        offsets[key] = read_long(fh_in)

    return offsets


with open("207513420127_R08C01_Grn.idat", "rb") as fh_in:
    print("1.  Magic:        ["+get_magic(fh_in)+"]")
    print("2.  IDAT version: ["+str(get_idat_version(fh_in))+"]")
    print("3.  section offsets:")
    
    offsets = get_section_offsets(fh_in)
    for key in sorted(offsets):
        if key in sections:
            print("    " + str(key) + ": " + str(offsets[key]) + "    ("+sections[key]+")")
        else:
            print("    " + str(key) + ": " + str(offsets[key]))

