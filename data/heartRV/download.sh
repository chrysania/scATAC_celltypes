cd ~/scratch/scATAC_celltypes/data/heartRV

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF223DBM/@@download/ENCFF223DBM.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF145ABZ/@@download/ENCFF145ABZ.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF400OJD/@@download/ENCFF400OJD.tar.gz #3
wget https://www.encodeproject.org/files/ENCFF075ZCA/@@download/ENCFF075ZCA.tar.gz #4
wget https://www.encodeproject.org/files/ENCFF358RDK/@@download/ENCFF358RDK.tar.gz #5
wget https://www.encodeproject.org/files/ENCFF962TJU/@@download/ENCFF962TJU.tar.gz #6
wget https://www.encodeproject.org/files/ENCFF781VVD/@@download/ENCFF781VVD.tar.gz #7
wget https://www.encodeproject.org/files/ENCFF717SDD/@@download/ENCFF717SDD.tar.gz #8
wget https://www.encodeproject.org/files/ENCFF655IKC/@@download/ENCFF655IKC.tar.gz #9
wget https://www.encodeproject.org/files/ENCFF537KXN/@@download/ENCFF537KXN.tar.gz #10
tar -xzvf ENCFF223DBM.tar.gz
tar -xzvf ENCFF145ABZ.tar.gz
tar -xzvf ENCFF400OJD.tar.gz
tar -xzvf ENCFF075ZCA.tar.gz
tar -xzvf ENCFF358RDK.tar.gz
tar -xzvf ENCFF962TJU.tar.gz
tar -xzvf ENCFF781VVD.tar.gz
tar -xzvf ENCFF717SDD.tar.gz
tar -xzvf ENCFF655IKC.tar.gz
tar -xzvf ENCFF537KXN.tar.gz

# RNA matrix (1-3 only)
mkdir -p heartRV1
cd heartRV1
wget https://www.encodeproject.org/files/ENCFF757DCI/@@download/ENCFF757DCI.tar.gz
tar -xzvf ENCFF757DCI.tar.gz
cd .. 

mkdir -p heartRV2
cd heartRV2
wget https://www.encodeproject.org/files/ENCFF289QCB/@@download/ENCFF289QCB.tar.gz
tar -xzvf ENCFF289QCB.tar.gz
cd .. 

mkdir -p heartRV3
cd heartRV2
wget https://www.encodeproject.org/files/ENCFF038SVX/@@download/ENCFF038SVX.tar.gz
tar -xzvf ENCFF038SVX.tar.gz
cd ..

mkdir -p heartRV4
mkdir -p heartRV5
mkdir -p heartRV6
mkdir -p heartRV7
mkdir -p heartRV8
mkdir -p heartRV9
mkdir -p heartRV10

# heartRV (right ventricle)
#sample_IDs=("ENCSR588PEE" "ENCSR681OLJ" "ENCSR169BCG" "ENCSR604PDO" "ENCSR814OLA" "ENCSR520ZUD" "ENCSR579TPC" "ENCSR615TSN" "ENCSR517QNQ" "ENCSR454YDZ")
#sample_n=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
#cell_type="heartRV"