idat-tools
==========

python toolkit to analyse, view and modify (mix) idat files.


Installation:

```{bash}
git clone https://github.com/yhoogstrate/idat-tools.git
virtualenv -p python3 .venv
source .venv/bin/activate
pip install .
idat-tools --version
```

## idat-tools view

Usage [idat-tools view]:
```{bash}
idat-tools view GSM6379997_203927450093_R01C01_Grn.idat 
```

Results in:

```{bash}
# manifest:             ''
# manifest (old style): ''
# unknown #1:           [1][0][0][0]
# sample id:            ''
# description:          ''
# plate:                ''
# well:                 ''
# unknown #2:           ''
# run info:
# 1. [8/21/2019 3:16:12 PM] [Decoding] [CallsToUsed=2798107|CallsToUnused=1889|CallsToInvalid=49292] [AutoDecode] [2.6.2]
# 2. [6/5/2020 2:44:03 PM] [Scan] [sherlockID=N1065|ScannerID=N1065|Username=Illumina|FPGAVersion=4.0.20|SoftwareAppication=iScan Control Software|SoftwareVersion=3.4.8] [iScan Control Software] [3.4.8]
# 3. [6/5/2020 2:44:03 PM] [Register] [Algorithm=StandardGeneric] [iScan Control Software] [3.4.8]
# 4. [6/5/2020 2:44:03 PM] [Extract] [Algorithm=StandardWithBackground] [iScan Control Software] [3.4.8]
# 5. [8/21/2019 3:16:26 PM] [Decoding] [CallsToUsed=2860469|CallsToUnused=1764|CallsToInvalid=60031] [AutoDecode] [2.6.2]
# 6. [6/5/2020 2:44:24 PM] [Scan] [sherlockID=N1065|ScannerID=N1065|Username=Illumina|FPGAVersion=4.0.20|SoftwareAppication=iScan Control Software|SoftwareVersion=3.4.8] [iScan Control Software] [3.4.8]
# 7. [6/5/2020 2:44:24 PM] [Register] [Algorithm=StandardGeneric] [iScan Control Software] [3.4.8]
# 8. [6/5/2020 2:44:24 PM] [Extract] [Algorithm=StandardWithBackground] [iScan Control Software] [3.4.8]
# 9. [8/21/2019 3:16:18 PM] [Decoding] [CallsToUsed=2867167|CallsToUnused=1847|CallsToInvalid=65541] [AutoDecode] [2.6.2]
# 10. [6/5/2020 2:44:44 PM] [Scan] [sherlockID=N1065|ScannerID=N1065|Username=Illumina|FPGAVersion=4.0.20|SoftwareAppication=iScan Control Software|SoftwareVersion=3.4.8] [iScan Control Software] [3.4.8]
# 11. [6/5/2020 2:44:44 PM] [Register] [Algorithm=StandardGeneric] [iScan Control Software] [3.4.8]
# 12. [6/5/2020 2:44:44 PM] [Extract] [Algorithm=StandardWithBackground] [iScan Control Software] [3.4.8]
# 13. [8/21/2019 3:16:10 PM] [Decoding] [CallsToUsed=2862808|CallsToUnused=1695|CallsToInvalid=53309] [AutoDecode] [2.6.2]
# 14. [6/5/2020 2:45:05 PM] [Scan] [sherlockID=N1065|ScannerID=N1065|Username=Illumina|FPGAVersion=4.0.20|SoftwareAppication=iScan Control Software|SoftwareVersion=3.4.8] [iScan Control Software] [3.4.8]
# 15. [6/5/2020 2:45:05 PM] [Register] [Algorithm=StandardGeneric] [iScan Control Software] [3.4.8]
# 16. [6/5/2020 2:45:05 PM] [Extract] [Algorithm=StandardWithBackground] [iScan Control Software] [3.4.8]
# 17. [8/21/2019 3:16:45 PM] [Decoding] [CallsToUsed=2886279|CallsToUnused=1866|CallsToInvalid=65578] [AutoDecode] [2.6.2]
# 18. [6/5/2020 2:45:24 PM] [Scan] [sherlockID=N1065|ScannerID=N1065|Username=Illumina|FPGAVersion=4.0.20|SoftwareAppication=iScan Control Software|SoftwareVersion=3.4.8] [iScan Control Software] [3.4.8]
# 19. [6/5/2020 2:45:24 PM] [Register] [Algorithm=StandardGeneric] [iScan Control Software] [3.4.8]
# 20. [6/5/2020 2:45:24 PM] [Extract] [Algorithm=StandardWithBackground] [iScan Control Software] [3.4.8]

IDAT v3: 203927450093_R01C01 (R/G: 0, BeadChip 8x5)
         probe_ids  probe_std_devs  probe_mean_intensities  probe_n_beads  probe_mid_block
0          1600101             436                    7582              9          1600101
1          1600111             523                    4757             14          1600111
2          1600115             554                    4851             18          1600115
3          1600123             299                    3165              9          1600123
4          1600131             152                     173             15          1600131
...            ...             ...                     ...            ...              ...
1051810   99810958             504                    2401             12         99810958
1051811   99810970              74                      97             10         99810970
1051812   99810978            1030                    7236             14         99810978
1051813   99810990             331                    4191              8         99810990
1051814   99810992             585                    4213             17         99810992

[1051815 rows x 5 columns]
```

## idat-tools view

Usage [idat-tools view]:
```{bash}
idat-tools mix --help
Usage: idat-tools mix [OPTIONS] IDAT_FILE_REFERENCE IDAT_FILE_MIXED_IN
                      IDAT_FILE_OUTPUT

Options:
  -r, --mix-ratio FLOAT RANGE  Fraction of mixed-in file values to be mixed
                               into reference file. E.g. 0.25 results in 75%
                               of reference and 25% of mixed-in file.
                               [default: 0.5; 0<=x<=1]
  --help                       Show this message and exit.
```
