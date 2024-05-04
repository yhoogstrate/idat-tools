#!/usr/bin/env python


from parser import *

for i in range(1000000):
    it1 = int_to_bytes(i, 8)
    it2 = bytes_to_int(it1, False) # set to true for long
    
    if i != it2:
        raise Exception("Err: "+str(i))


import sys
sys.exit(1)

# bytes_to_int(fh_in.read(8), signed=True)
# int_to_bytes


with open("GSM6379997_203927450093_R01C01_Grn.idat", "rb") as fh_1:
    with open("test.idat", "rb") as fh_2:
        for i in range(13676211):
            c1 = fh_1.read(1)
            c2 = fh_2.read(1)
            
            if c1 != c2:
                print("!! "+str(i) + "  " + str(c1) + " <> " + str(c2))
