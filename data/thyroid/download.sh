cd ~/scratch/scATAC_celltypes/data/thyroid

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF895JFF/@@download/ENCFF895JFF.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF891MZZ/@@download/ENCFF891MZZ.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF410LVV/@@download/ENCFF410LVV.tar.gz #3
tar -xzvf ENCFF895JFF.tar.gz
tar -xzvf ENCFF891MZZ.tar.gz
tar -xzvf ENCFF410LVV.tar.gz

mkdir -p thyroid1
mkdir -p thyroid2
mkdir -p thyroid3

# thyroid
#sample_IDs=("ENCSR796RXX-1" "ENCSR909OXO-1" "ENCSR817VFO-1")
#sample_n=("1" "2" "3")
#cell_type="thyroid"