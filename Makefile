

test:
	wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6379nnn/GSM6379997/suppl/GSM6379997%5F203927450093%5FR01C01%5FGrn.idat.gz ;
	wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6379nnn/GSM6379997/suppl/GSM6379997%5F203927450093%5FR01C01%5FRed.idat.gz ;
	wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5720495&format=file&file=GSM5720495%5F10003886252%5FR01C01%5FGrn%2Eidat%2Egz" ;
	wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5720495&format=file&file=GSM5720495%5F10003886252%5FR01C01%5FRed%2Eidat%2Egz" ;
	gunzip *.idat.gz

