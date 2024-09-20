cd ~/scratch/scATAC_celltypes/data/uterus

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF403BDC/@@download/ENCFF403BDC.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF215DOT/@@download/ENCFF215DOT.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF881QNT/@@download/ENCFF881QNT.tar.gz #3
tar -xzvf ENCFF403BDC.tar.gz
tar -xzvf ENCFF215DOT.tar.gz
tar -xzvf ENCFF881QNT.tar.gz

# RNA matrix (1 only)
mkdir -p uterus1
#wget https://www.encodeproject.org/files/ENCFF962ATV/@@download/ENCFF962ATV.tar.gz

mkdir -p uterus2
mkdir -p uterus3

# uterus
#sample_IDs=("ENCSR828HVB" "ENCSR455CVZ" "ENCSR028RFK")
#sample_n=("1" "2" "3")
#cell_type="uterus"