cd ~/scratch/scATAC_celltypes/data/leg_muscle

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

mkdir -p leg_muscle1
mkdir -p leg_muscle2
mkdir -p leg_muscle3
mkdir -p leg_muscle4
mkdir -p leg_muscle5

# gastrocnemius medialis (leg_muscle)
#sample_IDs=("ENCSR696YOC-1" "ENCSR819EGE-1" "ENCSR244GZL-1" "ENCSR023FME-1" "ENCSR139TIQ-1")
#sample_n=("1" "2" "3" "4" "5")
#cell_type="leg_muscle"