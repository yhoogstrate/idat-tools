#!/bin/bash python

with open("GSM6379997_203927450093_R01C01_Grn.idat", "rb") as fh_1:
    with open("test.idat", "rb") as fh_2:
        for i in range(13676211):
            c1 = fh_1.read(1)
            c2 = fh_2.read(1)
            
            if c1 != c2:
                print("!! "+str(i))
