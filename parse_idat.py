#!/usr/bin/env python

# https://github.com/FoxoTech/methylprep/blob/master/methylprep/files/idat.py#L228
# https://github.com/HenrikBengtsson/illuminaio/blob/develop/R/readIDAT_nonenc.R
# https://github.com/bioinformed/glu-genetics/blob/dcbbbf67a308d35e157b20a9c76373530510379a/glu/lib/illumina.py#L44-L61


from parser import *
from pathlib import Path
import re

from beartype import beartype
from _io import BufferedReader
from _io import BufferedWriter

import numpy as np
from numpy import ndarray

import pandas as pd
from pandas import DataFrame



section_names = {
    102: 'PROBE_IDS',
    103: 'PROBE_STD_DEVS',
    104: 'PROBE_MEAN_INTENSITIES',
    107: 'PROBE_N_BEADS',
    200: 'PROBE_MID_BLOCK', # also contains the probe_id's for some reason?
    300: 'ARRAY_RUN_INFO',
    400: 'ARRAY_RED_GREEN', # really concerned about this one, always [0]
    401: 'ARRAY_MANIFEST',
    402: 'ARRAY_BARCODE', # e.g. '203927450093'
    403: 'ARRAY_CHIP_TYPE', # e.g. 'BeadChip 8x5'
    404: 'ARRAY_CHIP_LABEL', # e.g. 'R01C01'
    405: 'ARRAY_OLD_STYLE_MANIFEST',
    406: 'ARRAY_SAMPLE_ID',
    407: 'ARRAY_DESCRIPTION',
    408: 'ARRAY_PLATE',
    409: 'ARRAY_WELL',
    410: 'ARRAY_UNKNOWN_1', # [1][0][0][0] <- could be int for 1, could be byte set, tuple of bytes is safest and easier to track new values down
    510: 'ARRAY_UNKNOWN_2',
    1000: 'ARRAY_N_PROBES'
}



class IDATdata(object):

    def __init__(self):
        self.data_file_magic = None
        self.data_idat_version = None
        self.data_section_order = None
        self.data_array_n_probes = None
        
        self.data_per_probe_matrix = None
        
        self.data_array_red_green = None
        self.data_array_manifest = None
        self.data_array_barcode = None
        self.data_array_chip_type = None
        self.data_array_chip_label = None
        self.data_array_old_style_manifest = None
        self.data_array_unknown_1 = None
        self.data_array_sample_id = None
        self.data_array_description = None
        self.data_array_plate = None
        self.data_array_well = None
        self.data_array_unknown_2 = None
        self.data_array_run_info: None
    
    def __str__(self):
        out = ""
        
        out += "# manifest:             '" + str(self.data_array_manifest) + "'\n"
        out += "# manifest (old style): '" + str(self.data_array_old_style_manifest) + "'\n"
        out += "# unknown #1:           [" + "][".join([str(_) for _ in self.data_array_unknown_1]) + "]\n"
        out += "# sample id:            '" + str(self.data_array_sample_id) + "'\n"
        out += "# description:          '" + str(self.data_array_description) + "'\n"
        out += "# plate:                '" + str(self.data_array_plate) + "'\n"
        out += "# well:                 '" + str(self.data_array_well) + "'\n"
        out += "# unknown #2:           '" + str(self.data_array_unknown_2) + "'\n"
        out += "# run info:\n"
        for i in range(len(self.data_array_run_info)):
            out += "# "+str(i+1)+". [" + "] [".join([str(_) for _ in self.data_array_run_info[i]]) + "]\n"

        out += "\n"
        out += self.data_file_magic 
        out += " v"
        out += str(self.data_idat_version)
        out += ": "
        out += self.data_array_barcode
        out += "_"
        out += self.data_array_chip_label
        out += " (R/G: "
        out += str(self.data_array_red_green)
        out += ", "
        out += self.data_array_chip_type
        out += ")"
        out += "\n"
        out += str(self.data_per_probe_matrix)

        return out
    
    @beartype
    def set_file_magic(self, file_magic: str) -> str:
        if file_magic != "IDAT":
            raise Exception("Invalid file format")
        else:
            self.data_file_magic = file_magic
        
        return self.data_file_magic
    
    @beartype
    def set_idat_version(self, idat_version: int) -> int:
        if idat_version != 3:
            raise Exception("Only tested with idat version 3")
        else:
            self.data_idat_version = idat_version
        
        return self.data_idat_version
    
    @beartype
    def set_array_n_probes(self, array_n_probes: int) -> int:
        if array_n_probes <= 0:
            raise Exception("Invalid number of probes: " + str(array_n_probes))
        else:
            self.data_array_n_probes = array_n_probes
        
        return self.data_array_n_probes

    @beartype
    def set_section_order(self, section_order: list[str]) -> list[str]:
        for _ in section_order:
            if _ not in section_names.values():
                raise Exception("Unknown section: "+str(_))
        
        self.data_section_order = section_order
        return self.data_section_order

    @beartype
    def set_data_per_probe_matrix(self, per_probe_matrix: DataFrame) -> DataFrame:
        if per_probe_matrix.shape[0] != self.data_array_n_probes:
            raise Exception("Matrix (nrow: "+str(per_probe_matrix.shape[0])+") does no fit size of the array (n="+str(self.data_array_n_probes)+")")
    
        if not per_probe_matrix['probe_ids'].equals(per_probe_matrix['probe_mid_block']):
            raise Exception("Discrepance between probe_ids and probe_mid_block")

        self.data_per_probe_matrix = per_probe_matrix
        return self.data_per_probe_matrix

    @beartype
    def set_array_red_green(self, array_red_green: int) -> int:
        # checks here
        
        self.data_array_red_green = array_red_green
        return self.data_array_red_green

    @beartype
    def set_array_manifest(self, array_manifest: str) -> str:
        # checks here
        
        self.data_array_manifest = array_manifest
        return self.data_array_manifest
    
    @beartype
    def set_array_barcode(self, array_barcode: str) -> str:
        if not re.match(r"^[0-9]+$", array_barcode):
            raise Exception("Incorrect barcode: " + array_barcode)
        
        self.data_array_barcode = array_barcode
        return self.data_array_barcode
    
    @beartype
    def set_array_chip_type(self, array_chip_type: str) -> str:
        if array_chip_type != "BeadChip 8x5":
            raise Exception("This tool is not tested with other chip types than BeadChip 8x5")
        
        self.data_array_chip_type = array_chip_type
        return self.data_array_chip_type
    
    @beartype
    def set_array_chip_label(self, array_chip_label: str) -> str:
        if not re.match(r"^R[0-9]+C[0-9]+$", array_chip_label):
            raise Exception("Odd label: " + array_chip_label)
        
        self.data_array_chip_label = array_chip_label
        return self.data_array_chip_label
    
    @beartype
    def set_array_old_style_manifest(self, array_old_style_manifest: str) -> str:
        # checks here
        
        self.data_array_old_style_manifest = array_old_style_manifest
        return self.data_array_old_style_manifest
    
    @beartype
    def set_array_unknown_1(self, array_unknown_1: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
        # checks here
        
        self.data_array_unknown_1 = array_unknown_1
        return self.data_array_unknown_1
    
    @beartype
    def set_array_sample_id(self, array_sample_id: str) -> str:
        # checks here
        
        self.data_array_sample_id = array_sample_id
        return self.data_array_sample_id

    @beartype
    def set_array_description(self, array_description: str) -> str:
        # checks here
        
        self.data_array_description = array_description
        return self.data_array_description

    @beartype
    def set_array_plate(self, array_plate: str) -> str:
        # checks here
        
        self.data_array_plate = array_plate
        return self.data_array_plate

    @beartype
    def set_array_well(self, array_well: str) -> str:
        # checks here
        
        self.data_array_well = array_well
        return self.data_array_well

    @beartype
    def set_array_unknown_2(self, array_unknown_2: str) -> str:
        # checks here
        
        self.data_array_unknown_2 = array_unknown_2
        return self.data_array_unknown_2
    
    @beartype
    def set_array_run_info(self, array_run_info: list[tuple[str, str, str, str, str]]) -> list[tuple[str, str, str, str, str]]:
        # checks here
        
        self.data_array_run_info = array_run_info
        return self.data_array_run_info



class IDATreader:

    @beartype
    def __init__(self, idat_filename: Path):
        self.data = IDATdata()
        
        self.idat_filename = idat_filename
        self.parse()
    
    @beartype
    def parse_file_magic(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['FILE_MAGIC'])
        
        return self.data.set_file_magic(read_char(fh_in, 4))
    
    @beartype
    def parse_idat_version(self, fh_in: BufferedReader, section_seek_index: dict) -> int:
        fh_in.seek(section_seek_index['IDAT_VERSION'])
        
        return self.data.set_idat_version(read_long(fh_in))
        
    @beartype
    def parse_section_index(self, fh_in: BufferedReader, section_seek_index: dict) -> dict:
        section_order = []
        
        fh_in.seek(section_seek_index['SECTION_INDEX_N'])
        n_sections = read_int(fh_in)
        
        for i in range(n_sections):
            section_type_int = read_short(fh_in)
            
            if section_type_int not in section_names:
                raise Exception("Unimplemented section type: "+str(section_type_int))
            else:
                section_type = section_names[section_type_int]
            
            section_order.append(section_type)
            section_seek_index[section_type] = read_long(fh_in)

        self.data.set_section_order(section_order)
        return section_seek_index
    
    @beartype
    def parse_array_n_probes(self, fh_in: BufferedReader, section_seek_index: dict) -> int:
        fh_in.seek(section_seek_index['ARRAY_N_PROBES'])
        
        return self.data.set_array_n_probes(read_int(fh_in))


    @beartype
    def parse_probe_ids(self, fh_in: BufferedReader, section_seek_index: dict) -> ndarray:
        fh_in.seek(section_seek_index['PROBE_IDS'])
        
        if self.data.data_array_n_probes is None:
            self.parse_array_n_probes(fh_in, section_seek_index)
        
        probe_ids = npread(fh_in, '<u4', self.data.data_array_n_probes) # layout-check: (4207470- 210) / 1051815 = 4
        
        if np.any(probe_ids <= 0):
            raise Exception("Wrong probe id's found")
        
        if np.any((probe_ids[1:] - probe_ids[:-1]) <= 0):
            raise Exception("probe id's are not unique or not incremental")
        
        return probe_ids

    @beartype
    def parse_probe_std_devs(self, fh_in: BufferedReader, section_seek_index: dict) -> ndarray:
        fh_in.seek(section_seek_index['PROBE_STD_DEVS'])
        
        if self.data.data_array_n_probes is None:
            self.parse_array_n_probes(fh_in, section_seek_index)
        
        probe_std_devs = npread(fh_in, '<u2', self.data.data_array_n_probes) # layout-check: (6311100 - 4207470) / 1051815 = 2
        
        if np.any(probe_std_devs < 0):
            raise Exception("Wrong std dev found (0 or negative)")
        
        return probe_std_devs

    @beartype
    def parse_probe_mean_intensities(self, fh_in: BufferedReader, section_seek_index: dict) -> ndarray:
        fh_in.seek(section_seek_index['PROBE_MEAN_INTENSITIES'])
        
        if self.data.data_array_n_probes is None:
            self.parse_array_n_probes(fh_in, section_seek_index)
        
        probe_mean_intensities = npread(fh_in, '<u2', self.data.data_array_n_probes) # layout-check: (8414730 - 6311100) / 1051815 = 2
        
        if np.any(probe_mean_intensities < 0):
            raise Exception("Wrong median probe intensity found (negative)")
        
        return probe_mean_intensities

    @beartype
    def parse_probe_n_beads(self, fh_in: BufferedReader, section_seek_index: dict) -> ndarray:
        fh_in.seek(section_seek_index['PROBE_N_BEADS'])
        
        if self.data.data_array_n_probes is None:
            self.parse_array_n_probes(fh_in, section_seek_index)
        
        probe_n_beads = npread(fh_in, '<u1', self.data.data_array_n_probes) # layout-check: (9466545 - 8414730) / 1051815 = 1
        
        if np.any(probe_n_beads < 0):
            raise Exception("Wrong number of beads per probe found (0 or negative)")
        
        return probe_n_beads

    @beartype
    def parse_probe_mid_block(self, fh_in: BufferedReader, section_seek_index: dict) -> ndarray:
        fh_in.seek(section_seek_index['PROBE_MID_BLOCK'])
        
        if self.data.data_array_n_probes is None:
            self.parse_array_n_probes(fh_in, section_seek_index)
        
        if self.data.data_array_n_probes != read_int(fh_in):
            raise Exception("Weird discrepancy between number of probes and size of mid block")
        
        probe_mid_block = npread(fh_in, '<u4', self.data.data_array_n_probes) # layout-check: (13673809 - (9466545 + 4)) / 1051815 = 4
        
        if np.any(probe_mid_block <= 0):
            raise Exception("Wrong probe id's found")
        
        if np.any((probe_mid_block[1:] - probe_mid_block[:-1]) <= 0):
            raise Exception("probe id's are not unique or not incremental")
        
        return probe_mid_block

    @beartype
    def parse_data_per_probe_matrix(self, fh_in: BufferedReader, section_seek_index: dict) -> DataFrame:
        data_per_probe_matrix = pd.DataFrame({
            'probe_ids': self.parse_probe_ids(fh_in, section_seek_index),
            'probe_std_devs': self.parse_probe_std_devs(fh_in, section_seek_index),
            'probe_mean_intensities': self.parse_probe_mean_intensities(fh_in, section_seek_index),
            'probe_n_beads': self.parse_probe_n_beads(fh_in, section_seek_index),
            'probe_mid_block': self.parse_probe_mid_block(fh_in, section_seek_index)
            })

        if not data_per_probe_matrix['probe_ids'].equals(data_per_probe_matrix['probe_mid_block']):
            raise Exception("Discrepance between probe_ids and probe_mid_block")

        self.data_per_probe_matrix = data_per_probe_matrix

        return self.data.set_data_per_probe_matrix(data_per_probe_matrix)

    @beartype
    def parse_array_red_green(self, fh_in: BufferedReader, section_seek_index: dict) -> int:
        fh_in.seek(section_seek_index['ARRAY_RED_GREEN'])
        
        red_green = read_int(fh_in)
        
        #if red_green != 0:
        #    raise Exception("Only seen 0 so far, but probably good...")
        
        return self.data.set_array_red_green(red_green)

    @beartype
    def parse_array_manifest(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_MANIFEST'])
        
        return self.data.set_array_manifest(read_string(fh_in))
    
    @beartype
    def parse_array_barcode(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_BARCODE'])
        
        return self.data.set_array_barcode(read_string(fh_in))

    @beartype
    def parse_array_chip_type(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_CHIP_TYPE'])
        
        return self.data.set_array_chip_type(read_string(fh_in))

    @beartype
    def parse_array_chip_label(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_CHIP_LABEL'])
        
        return self.data.set_array_chip_label(read_string(fh_in))

    @beartype
    def parse_array_old_style_manifest(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_OLD_STYLE_MANIFEST'])
        
        return self.data.set_array_old_style_manifest(read_string(fh_in))

    @beartype
    def parse_array_unknown_1(self, fh_in: BufferedReader, section_seek_index: dict) -> tuple[int, int, int, int]:
        fh_in.seek(section_seek_index['ARRAY_UNKNOWN_1'])
        
        return self.data.set_array_unknown_1((read_byte(fh_in), read_byte(fh_in), read_byte(fh_in), read_byte(fh_in)))

    @beartype
    def parse_array_sample_id(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_SAMPLE_ID'])
        
        return self.data.set_array_sample_id(read_string(fh_in))

    @beartype
    def parse_array_description(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_DESCRIPTION'])
        
        return self.data.set_array_description(read_string(fh_in))

    @beartype
    def parse_array_plate(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_PLATE'])
        
        return self.data.set_array_plate(read_string(fh_in))

    @beartype
    def parse_array_well(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_WELL'])
        
        return self.data.set_array_well(read_string(fh_in))

    @beartype
    def parse_array_unknown_2(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_UNKNOWN_2'])
        
        return self.data.set_array_unknown_2(read_string(fh_in))
    
    @beartype
    def parse_array_run_info(self, fh_in: BufferedReader, section_seek_index: dict) -> list[tuple[str, str, str, str, str]]:
        fh_in.seek(section_seek_index['ARRAY_RUN_INFO'])
        
        run_info = []
        
        for i in range(read_int(fh_in)): # blocks containing 5 consequtive strings
            run_info.append((read_string(fh_in), read_string(fh_in), read_string(fh_in), read_string(fh_in), read_string(fh_in)))
        
        return self.data.set_array_run_info(run_info)

    @beartype
    def parse(self) -> int:
        section_seek_index = { # static entries - todo: make class
            'FILE_MAGIC': 0,
            'IDAT_VERSION': 4,
            'SECTION_INDEX_N': 12
        }
        
        with open(self.idat_filename, "rb") as fh_in:
            self.parse_file_magic(fh_in, section_seek_index)
            self.parse_idat_version(fh_in, section_seek_index)
            self.parse_section_index(fh_in, section_seek_index)
            self.parse_array_n_probes(fh_in, section_seek_index)
            self.parse_data_per_probe_matrix(fh_in, section_seek_index)
            self.parse_array_red_green(fh_in, section_seek_index)
            self.parse_array_manifest(fh_in, section_seek_index)
            self.parse_array_barcode(fh_in, section_seek_index)
            self.parse_array_chip_type(fh_in, section_seek_index)
            self.parse_array_chip_label(fh_in, section_seek_index)
            self.parse_array_old_style_manifest(fh_in, section_seek_index)
            self.parse_array_unknown_1(fh_in, section_seek_index)
            self.parse_array_sample_id(fh_in, section_seek_index)
            self.parse_array_description(fh_in, section_seek_index)
            self.parse_array_plate(fh_in, section_seek_index)
            self.parse_array_well(fh_in, section_seek_index)
            self.parse_array_unknown_2(fh_in, section_seek_index)
            self.parse_array_run_info(fh_in, section_seek_index)

        return 0



class IDATwriter(IDATdata):

    @beartype
    def __init__(self, idat_data: IDATdata):
        if isinstance(idat_data, IDATdata):
            self.data = idat_data
        elif isinstance(idat_data, IDATreader):
            self.data = idat_data.data
        elif isinstance(idat_data, IDATwriter):
            self.data = idat_data.data
        else:
            raise Exception("Unclear input type")

    @beartype
    def write(self, idat_filename: Path):
        with open(idat_filename, 'wb') as fh_out:
            write_char(fh_out, self.data.data_file_magic)
            write_long(fh_out, self.data.data_idat_version)



d_red = IDATreader(Path("GSM6379997_203927450093_R01C01_Grn.idat"))

w = IDATwriter(d_red.data)
w.write(Path("test.idat"))

#d_grn = IDATreader(Path("GSM6379997_203927450093_R01C01_Red.idat"))




