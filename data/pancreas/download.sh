cd ~/scratch/scATAC_celltypes/data/pancreas

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF662CFN/@@download/ENCFF662CFN.tar.gz
wget https://www.encodeproject.org/files/ENCFF836KOE/@@download/ENCFF836KOE.tar.gz
wget https://www.encodeproject.org/files/ENCFF482UTI/@@download/ENCFF482UTI.tar.gz
wget https://www.encodeproject.org/files/ENCFF455WME/@@download/ENCFF455WME.tar.gz
tar -xzvf ENCFF662CFN.tar.gz
tar -xzvf ENCFF836KOE.tar.gz
tar -xzvf ENCFF482UTI.tar.gz
tar -xzvf ENCFF455WME.tar.gz

# RNA matrix
mkdir -p pancreas1
cd pancreas1
wget https://www.encodeproject.org/files/ENCFF083VEZ/@@download/ENCFF083VEZ.tar.gz
tar -xzvf ENCFF083VEZ.tar.gz
cd ..

mkdir pancreas2
cd pancreas2
wget https://www.encodeproject.org/files/ENCFF069UYQ/@@download/ENCFF069UYQ.tar.gz
tar -xzvf ENCFF069UYQ.tar.gz
cd ..

mkdir pancreas3
cd pancreas3
wget https://www.encodeproject.org/files/ENCFF768ULN/@@download/ENCFF768ULN.tar.gz
tar -xzvf ENCFF768ULN.tar.gz
cd ..

mkdir pancreas4
cd pancreas4
wget https://www.encodeproject.org/files/ENCFF012YZW/@@download/ENCFF012YZW.tar.gz
tar -xzvf ENCFF012YZW.tar.gz
cd ..

# pancreas
#sample_IDs=("ENCSR868CRK-1" "ENCSR690NZI-1" "ENCSR229VVY-1" "ENCSR496XXB-1")
#sample_n=("1" "2" "3" "4")
#cell_type="pancreas"