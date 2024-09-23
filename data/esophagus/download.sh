cd ~/scratch/scATAC_celltypes/data/esophagus

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF026ZEV/@@download/ENCFF026ZEV.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF623PSE/@@download/ENCFF623PSE.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF815ESF/@@download/ENCFF815ESF.tar.gz #3
tar -xzvf ENCFF026ZEV.tar.gz
tar -xzvf ENCFF623PSE.tar.gz
tar -xzvf ENCFF815ESF.tar.gz

mkdir -p esophagus1
mkdir -p esophagus2
mkdir -p esophagus3

# esophagus
#sample_IDs=("ENCSR453TVZ-1" "ENCSR757EGB-1" "ENCSR164GSH-1")
#sample_n=("1" "2" "3")
#cell_type="esophagus"