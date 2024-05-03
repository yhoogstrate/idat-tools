#!/usr/bin/env python

import numpy as np

def bytes_to_int(input_bytes, signed=False):
    """Returns the integer represented by the given array of bytes.
    Pre-sets the byteorder to be little-endian.

    Arguments:
        input_bytes -- Holds the array of bytes to convert.  The argument must either
            support the buffer protocol or be an iterable object producing bytes.
            Bytes and bytearray are examples of built-in objects that support the
            buffer protocol.

    Keyword Arguments:
        signed {bool} -- Indicates whether two's complement is used to represent the integer. (default: {False})

    Returns:
        [integer] -- Integer value converted from the supplied bytes.
    """
    return int.from_bytes(input_bytes, byteorder='little', signed=signed)



def read_byte(infile):
    """Converts a single byte to an integer value.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Unsigned integer value converted from the supplied bytes.
    """
    return bytes_to_int(infile.read(1), signed=False)


def read_short(infile):
    """Converts a two-byte element to an integer value.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Unsigned integer value converted from the supplied bytes.
    """
    return bytes_to_int(infile.read(2), signed=False)


def read_int(infile):
    """Converts a four-byte element to an integer value.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Signed integer value converted from the supplied bytes.
    """
    return bytes_to_int(infile.read(4), signed=True)


def read_long(infile):
    """Converts an eight-byte element to an integer value.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Signed integer value converted from the supplied bytes.
    """
    return bytes_to_int(infile.read(8), signed=True)


def read_char(infile, num_bytes):
    """Converts an array of bytes to a string.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.
        num_bytes {integer} -- The number of bytes to read and parse.

    Returns:
        [string] -- UTF-8 decoded string value.
    """
    return infile.read(num_bytes).decode('utf-8')


def read_string(infile):
    """Converts an array of bytes to a string.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [string] -- UTF-8 decoded string value.
    """
    num_bytes = read_byte(infile)
    num_chars = num_bytes % 128
    shift = 0

    while num_bytes // 128 == 1:
        num_bytes = read_byte(infile)
        shift += 7
        offset = (num_bytes % 128) * (2 ** shift)
        num_chars += offset

    return read_char(infile, num_chars)

def npread(file_like, dtype, n):
    # https://stackoverflow.com/questions/72838939/how-to-convert-the-string-between-numpy-array-and-bytes
    """Parses a binary file multiple times, allowing for control if the
    file ends prematurely. This replaces read_results() and runs faster.
    And it provides support for reading gzipped idat files without decompressing.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.
        dtype {data type} -- used within idat files, 2-bit, or 4-bit numbers stored in binary at specific addresses
        n {number of snps read} -- see files/idat.py for how this function is applied.

    Raises:
        EOFError: If the end of the file is reached before the number of elements have
            been processed.

    Returns:
        A list of the parsed values.
    """
    dtype = np.dtype(dtype)
    
    # np.readfile is not able to read from gzopene-d file
    alldata = file_like.read(dtype.itemsize * n)
    
    if len(alldata) != dtype.itemsize * n:
        raise EOFError('End of file reached before number of results parsed')
    readdata=np.frombuffer(alldata, dtype, n)
    if readdata.size != n:
        raise EOFError('End of file reached before number of results parsed')
    
    """
    # cast back and check
    result_str = np.ndarray.tobytes(readdata) #.decode("utf-8")
    
    print("alldata", alldata[1:5])
    print("results", result_str[1:5])
    if alldata[1:5] == result_str[1:5]:
        print("check , same!")
    """
    
    return readdata


