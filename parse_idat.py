#!/usr/bin/env python

# https://github.com/FoxoTech/methylprep/blob/master/methylprep/files/idat.py#L228
# https://github.com/HenrikBengtsson/illuminaio/blob/develop/R/readIDAT_nonenc.R
# https://github.com/bioinformed/glu-genetics/blob/dcbbbf67a308d35e157b20a9c76373530510379a/glu/lib/illumina.py#L44-L61


from parser import *
from beartype import beartype
import pandas as pd
from pathlib import Path
from _io import BufferedReader

section_names = {
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

section_locations = {
    'ILLUMINA_ID': 102,
    'STD_DEV': 103,
    'MEAN': 104,
    'NUM_BEADS': 107,
    'MID_BLOCK': 200,
    'RUN_INFO': 300,
    'RED_GREEN': 400,
    'MOSTLY_NULL': 401,
    'BARCODE': 402,
    'CHIP_TYPE': 403,
    'MOSTLY_A': 404,
    'UNKNOWN_1': 405,
    'UNKNOWN_2': 406,
    'UNKNOWN_3': 407,
    'UNKNOWN_4': 408,
    'UNKNOWN_5': 409,
    'UNKNOWN_6': 410,
    'UNKNOWN_7': 510,
    'NUM_SNPS_READ': 1000
}




class IDATdata:
    def __init__(self):
        self.data_version = None
        self.data_section_order = None
        self.data_n_probes = None
        self.data_sections = {
            'RED_GREEN': None,
            'MOSTLY_NULL': None,
            'BARCODE': None,
            'CHIP_TYPE': None,
            'MOSTLY_A': None,
            'UNKNOWN_1': None,
            'UNKNOWN_6': None,
            'UNKNOWN_2': None,
            'UNKNOWN_3': None,
            'UNKNOWN_4': None,
            'UNKNOWN_5': None,
            'UNKNOWN_7': None,
            'RUN_INFO': None
        }



class IDATfile(IDATdata):

    @beartype
    def __init__(self, idat_filename: Path):
        self.idat_filename = idat_filename
        self.parse()
    
    @beartype
    def parse_magic(self, fh_in: BufferedReader, section_seek_index):
        print(type(fh_in))
    
    @beartype
    def parse(self) -> int:
        section_seek_index = {
            'FILE_MAGIC': 0,
            'IDAT_VERSION': 4,
            'SECTION_INDEX_N': 12,
            'SECTION_INDEX': 16
        }
        
        with open(self.idat_filename, "rb") as fh_in:
            self.parse_magic(fh_in, section_seek_index)
        
        return 0


d_red = IDATfile(Path("207513420108_R01C01_Red.idat"))
d_grn = IDATfile(Path("207513420108_R01C01_Grn.idat"))





def get_magic(fh_in):
    fh_in.seek(0) # IdatHeaderLocation.FILE_TYPE.value
    file_type = read_char(fh_in, 4)
    
    return file_type

def get_idat_version(fh_in):
    fh_in.seek(4) # IdatHeaderLocation.VERSION.value
    print("    => [" + str(read_byte(fh_in))+"][" + str(read_byte(fh_in))+"][" + str(read_byte(fh_in))+"][" + str(read_byte(fh_in))+"][" + str(read_byte(fh_in))+"][" + str(read_byte(fh_in))+"][" + str(read_byte(fh_in))+"][" + str(read_byte(fh_in))+"]")

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
    pos = 16
    poss = []

    offsets = {}
    for _idx in range(num_fields):
        poss.append(str(pos) + "-" + str(pos+9))
        
        key = read_short(fh_in)
        offsets[read_long(fh_in)] = key
        
        pos += 2 + 8

    print("    => "+str(poss))

    return offsets


with open(Path("207513420108_R01C01_Red.idat"), "rb") as fh_in:
    print("1.  Magic:        ["+get_magic(fh_in)+"]")
    print("2.  IDAT version: ["+str(get_idat_version(fh_in))+"]")
    print("3.  section offsets:")
    
    offsets = get_section_offsets(fh_in)
    
    print("4.  details")
    for offset in sorted(offsets):
        key = offsets[offset]
        
        if key in section_names:
            print("    " + str(key) + ": " + str(offset) + "    ("+section_names[key]+")")
        else:
            print("    " + str(key) + ": " + str(offset))
        
        if key == 102:
            pass # just value parsing
        elif key == 103:
            pass # just value parsing
        elif key == 104:
            pass # just value parsing
        elif key == 107:
            pass # just value parsing
            #seek_to_section(IdatSectionCode.NUM_BEADS)
            #self.n_beads = npread(idat_file, '<u1', self.n_snps_read) # was <u1

        elif key == 200:
            fh_in.seek(offset)
            midblock = []
            for i in range(read_int(fh_in)):
                midblock.append(read_int(fh_in))
            
            midblock = np.array(midblock)
            print("    => len: " + str(len(midblock)) + ",  [" + str(midblock[0]) + ", "+str(midblock[1]) + ", ..., "  + str(midblock[-2]) + ", "+str(midblock[-1]) + "]")


        elif key == 300:
            fh_in.seek(offset)
            for i in range(read_int(fh_in)):
                print("    => "+str(i+1)+". [" + str(read_string(fh_in))+"][" + str(read_string(fh_in))+"][" + str(read_string(fh_in))+"][" + str(read_string(fh_in))+"][" + str(read_string(fh_in))+"]")
            
        elif key == 400:
            fh_in.seek(offset)
            print("    => [" + str(read_int(fh_in))+"]")
        elif key == 401:
            fh_in.seek(offset)
            print("    => [" + str(read_byte(fh_in))+"]")
        elif key == 402:
            fh_in.seek(offset)
            print("    => [" + str(read_string(fh_in))+"]")
        elif key == 403:
            fh_in.seek(offset)
            print("    => [" + str(read_string(fh_in))+"]")
        elif key == 404:
            fh_in.seek(offset)
            print("    => [" + str(read_string(fh_in))+"]")
        elif key == 405:
            fh_in.seek(offset)
            print("    => [" + str(read_byte(fh_in))+"]")
        elif key == 406:
            fh_in.seek(offset)
            print("    => [" + str(read_byte(fh_in))+"]")
        elif key == 407:
            fh_in.seek(offset)
            print("    => [" + str(read_byte(fh_in))+"]")
        elif key == 408:
            fh_in.seek(offset)
            print("    => [" + str(read_byte(fh_in))+"]")
        elif key == 409:
            fh_in.seek(offset)
            print("    => [" + str(read_byte(fh_in))+"]")
        elif key == 410:
            fh_in.seek(offset)
            print("    => [" + str(read_byte(fh_in))+"][" + str(read_byte(fh_in))+"][" + str(read_byte(fh_in))+"][" + str(read_byte(fh_in))+"]")
        elif key == 510:
            fh_in.seek(offset)
            print("    => [" + str(read_byte(fh_in))+"]")
        elif key == 1000:
            fh_in.seek(offset)
            print("    => [" + str(read_int(fh_in))+"]")
        else:
            raise Exception("Key: " + str(key) + " not yet implemented")


    offset = [_ for _ in offsets if offsets[_] == section_locations['NUM_SNPS_READ']][0]
    fh_in.seek(offset)
    n_snps_read = read_int(fh_in)

    offset = [_ for _ in offsets if offsets[_] == section_locations['NUM_BEADS']][0]
    fh_in.seek(offset)
    n_beads = npread(fh_in, '<u1', n_snps_read) # was <u1
    
    print("5. beads:  " + str(n_beads[0:24]) + "    (n="+str(len(n_beads))+")")
    
    
    offset = [_ for _ in offsets if offsets[_] == section_locations['ILLUMINA_ID']][0]
    fh_in.seek(offset)
    illumn = npread(fh_in, '<i4', n_snps_read) # was <u1
    
    print("6. Illumina IDs:  " + str(illumn[0:8]) + "    (n="+str(len(illumn))+")")
    print("   [midblock]  :  " + str(midblock[0:8]) + "    (n="+str(len(midblock))+")")
    
    if np.all(illumn == midblock):
        print("   => idenictal")
    else:
        raise Exception("discrepancies between identifier blocks")


    offset = [_ for _ in offsets if offsets[_] == section_locations['MEAN']][0]
    fh_in.seek(offset)
    mean = npread(fh_in, '<i4', n_snps_read) # was <u1
    
    print("7. Means:  " + str(mean[0:7]) + "    (n="+str(len(mean))+")")
    
    
    offset = [_ for _ in offsets if offsets[_] == section_locations['STD_DEV']][0]
    fh_in.seek(offset)
    sd = npread(fh_in, '<i4', n_snps_read) # was <u1
    
    print("8. Stddevs:  " + str(sd[0:7]) + "    (n="+str(len(sd))+")")
    
    
    df = pd.DataFrame({
                'probe_id': illumn,
                'probe_n_beads': n_beads,
                'probe_mean': mean,
                'probe_sd': sd })
    
    print(df)



