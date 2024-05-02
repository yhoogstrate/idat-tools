

test:
	wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6379nnn/GSM6379997/suppl/GSM6379997%5F203927450093%5FR01C01%5FGrn.idat.gz ;
	wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6379nnn/GSM6379997/suppl/GSM6379997%5F203927450093%5FR01C01%5FRed.idat.gz ;
	gunzip *.idat.gz
