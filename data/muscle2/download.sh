cd ~/scratch/scATAC_celltypes/data/muscle2

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF469KDC/@@download/ENCFF469KDC.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF738CZX/@@download/ENCFF738CZX.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF982BFF/@@download/ENCFF982BFF.tar.gz #3
wget https://www.encodeproject.org/files/ENCFF223JJZ/@@download/ENCFF223JJZ.tar.gz #4
wget https://www.encodeproject.org/files/ENCFF623EVG/@@download/ENCFF623EVG.tar.gz #5
tar -xzvf ENCFF469KDC.tar.gz
tar -xzvf ENCFF738CZX.tar.gz
tar -xzvf ENCFF982BFF.tar.gz
tar -xzvf ENCFF223JJZ.tar.gz
tar -xzvf ENCFF623EVG.tar.gz

mkdir -p muscle1
mkdir -p muscle2
mkdir -p muscle3
mkdir -p muscle4
mkdir -p muscle5

# gastrocnemius medialis (muscle2)
#sample_IDs=("ENCSR696YOC" "ENCSR819EGE" "ENCSR244GZL" "ENCSR023FME" "ENCSR139TIQ")
#sample_n=("1" "2" "3" "4" "5")
#cell_type="muscle2"