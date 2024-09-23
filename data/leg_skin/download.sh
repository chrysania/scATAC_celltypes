cd ~/scratch/scATAC_celltypes/data/leg_skin

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF327EKU/@@download/ENCFF327EKU.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF558THN/@@download/ENCFF558THN.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF299IOQ/@@download/ENCFF299IOQ.tar.gz #3
wget https://www.encodeproject.org/files/ENCFF785OMO/@@download/ENCFF785OMO.tar.gz #4
tar -xzvf ENCFF327EKU.tar.gz
tar -xzvf ENCFF558THN.tar.gz
tar -xzvf ENCFF299IOQ.tar.gz
tar -xzvf ENCFF785OMO.tar.gz

mkdir -p leg_skin1
mkdir -p leg_skin2
mkdir -p leg_skin3
mkdir -p leg_skin4

# lower leg skin
#sample_IDs=("ENCSR733SZL-1" "ENCSR397ODX-1" "ENCSR474TGL-1" "ENCSR513HZN-1")
#sample_n=("1" "2" "3" "4")
#cell_type="leg_skin"