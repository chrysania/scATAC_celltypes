cd ~/scratch/scATAC_celltypes/data/tibial_nerve

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF998LXW/@@download/ENCFF998LXW.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF161PNR/@@download/ENCFF161PNR.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF462EMP/@@download/ENCFF462EMP.tar.gz #3
tar -xzvf ENCFF998LXW.tar.gz
tar -xzvf ENCFF161PNR.tar.gz
tar -xzvf ENCFF462EMP.tar.gz

mkdir -p tibial_nerve1
mkdir -p tibial_nerve2
mkdir -p tibial_nerve3

# tibial nerve
#sample_IDs=("ENCSR205TUH-1" "ENCSR453BVR-1" "ENCSR726QTF-1")
#sample_n=("1" "2" "3")
#cell_type="tibial_nerve"