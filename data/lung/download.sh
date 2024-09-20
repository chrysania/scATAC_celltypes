cd ~/scratch/scATAC_celltypes/data/lung

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF184CZI/@@download/ENCFF184CZI.tar.gz
wget https://www.encodeproject.org/files/ENCFF175ZOV/@@download/ENCFF175ZOV.tar.gz
wget https://www.encodeproject.org/files/ENCFF280ZNE/@@download/ENCFF280ZNE.tar.gz
tar -xzvf ENCFF184CZI.tar.gz
tar -xzvf ENCFF175ZOV.tar.gz
tar -xzvf ENCFF280ZNE.tar.gz

# RNA matrix
mkdir -p lung1
cd lung1
wget https://www.encodeproject.org/files/ENCFF313RZL/@@download/ENCFF313RZL.tar.gz
tar -xzvf ENCFF313RZL.tar.gz
cd ..

mkdir -p lung2
cd lung2
wget https://www.encodeproject.org/files/ENCFF904BER/@@download/ENCFF904BER.tar.gz
tar -xzvf ENCFF904BER.tar.gz
cd ..

mkdir -p lung3
cd lung3
wget https://www.encodeproject.org/files/ENCFF415KNH/@@download/ENCFF415KNH.tar.gz
tar -xzvf ENCFF415KNH.tar.gz
cd ..

# lung
#sample_IDs=("ENCSR816NWE-1" "ENCSR391BWM-1" "ENCSR824OCY-1")
#sample_n=("1" "2" "3")
#cell_type="lung"