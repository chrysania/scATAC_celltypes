cd ~/scratch/scATAC_celltypes/data/left_colon

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF040VYY/@@download/ENCFF040VYY.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF810NQZ/@@download/ENCFF810NQZ.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF826JKT/@@download/ENCFF826JKT.tar.gz #3
tar -xzvf ENCFF040VYY.tar.gz
tar -xzvf ENCFF810NQZ.tar.gz
tar -xzvf ENCFF826JKT.tar.gz

# RNA matrix (1 only)
mkdir -p left_colon1
#https://www.encodeproject.org/files/ENCFF949IOU/@@download/ENCFF949IOU.tar.gz
mkdir -p left_colon2
mkdir -p left_colon3

# left colon
#sample_IDs=("ENCSR830FPR-1" "ENCSR916RYB-1" "ENCSR904WIW-1")
#sample_n=("1" "2" "3")
#cell_type="left_colon"