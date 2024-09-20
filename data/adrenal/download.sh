cd ~/scratch/scATAC_celltypes/data/adrenal

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF318DQR/@@download/ENCFF318DQR.tar.gz
wget https://www.encodeproject.org/files/ENCFF370JZO/@@download/ENCFF370JZO.tar.gz
wget https://www.encodeproject.org/files/ENCFF035BPL/@@download/ENCFF035BPL.tar.gz
tar -xzvf ENCFF318DQR.tar.gz
tar -xzvf ENCFF370JZO.tar.gz
tar -xzvf ENCFF035BPL.tar.gz

# RNA matrix
mkdir -p adrenal1
cd adrenal1
wget https://www.encodeproject.org/files/ENCFF451NNZ/@@download/ENCFF451NNZ.tar.gz
tar -xzvf ENCFF451NNZ.tar.gz
cd ..

mkdir -p adrenal2
cd adrenal2
wget https://www.encodeproject.org/files/ENCFF846BGR/@@download/ENCFF846BGR.tar.gz
tar -xzvf ENCFF846BGR.tar.gz
cd ..

mkdir -p adrenal3
cd adrenal3
wget https://www.encodeproject.org/files/ENCFF508AYA/@@download/ENCFF508AYA.tar.gz
tar -xzvf ENCFF508AYA.tar.gz
cd ..

# adrenal gland
#sample_IDs=("ENCSR420EWQ-1" "ENCSR693GAD-1" "ENCSR194KHA-1")
#sample_n=("1" "2" "3")
#cell_type="adrenal"