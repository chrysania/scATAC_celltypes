cd ~/scratch/scATAC_celltypes/data/omental_fat

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF261CAT/@@download/ENCFF261CAT.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF115LUT/@@download/ENCFF115LUT.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF768WMO/@@download/ENCFF768WMO.tar.gz #3
wget https://www.encodeproject.org/files/ENCFF614TFZ/@@download/ENCFF614TFZ.tar.gz #4
tar -xzvf ENCFF261CAT.tar.gz
tar -xzvf ENCFF115LUT.tar.gz
tar -xzvf ENCFF768WMO.tar.gz
tar -xzvf ENCFF614TFZ.tar.gz

mkdir -p omental_fat1
mkdir -p omental_fat2
mkdir -p omental_fat3
mkdir -p omental_fat4

# omental fat pad
#sample_IDs=("ENCSR181XXQ-1" "ENCSR492GGN-1" "ENCSR644SCP-1" "ENCSR274HQD-1")
#sample_n=("1" "2" "3" "4")
#cell_type="omental_fat"