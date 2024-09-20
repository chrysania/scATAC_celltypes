cd ~/scratch/scATAC_celltypes/data/liver

# ATAC matrix
#wget https://www.encodeproject.org/files/ENCFF478ULU/@@download/ENCFF478ULU.tar.gz #1
#wget https://www.encodeproject.org/files/ENCFF375OTA/@@download/ENCFF375OTA.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF589YDE/@@download/ENCFF589YDE.tar.gz #3
#wget https://www.encodeproject.org/files/ENCFF744GTE/@@download/ENCFF744GTE.tar.gz #4
#wget https://www.encodeproject.org/files/ENCFF033LOH/@@download/ENCFF033LOH.tar.gz #5
#wget https://www.encodeproject.org/files/ENCFF104HAE/@@download/ENCFF104HAE.tar.gz #6
wget https://www.encodeproject.org/files/ENCFF883DOK/@@download/ENCFF883DOK.tar.gz #7
#wget https://www.encodeproject.org/files/ENCFF397EDU/@@download/ENCFF397EDU.tar.gz #8
#wget https://www.encodeproject.org/files/ENCFF673SPK/@@download/ENCFF673SPK.tar.gz #9
wget https://www.encodeproject.org/files/ENCFF240WBA/@@download/ENCFF240WBA.tar.gz #10
#tar -xzvf ENCFF478ULU.tar.gz
#tar -xzvf ENCFF375OTA.tar.gz
tar -xzvf ENCFF589YDE.tar.gz
#tar -xzvf ENCFF744GTE.tar.gz
#tar -xzvf ENCFF033LOH.tar.gz
#tar -xzvf ENCFF104HAE.tar.gz
tar -xzvf ENCFF883DOK.tar.gz
#tar -xzvf ENCFF397EDU.tar.gz
#tar -xzvf ENCFF673SPK.tar.gz
tar -xzvf ENCFF240WBA.tar.gz

# RNA matrix
#mkdir -p liver1
#cd liver1
#wget https://www.encodeproject.org/files/ENCFF507QES/@@download/ENCFF507QES.tar.gz
#tar -xzvf ENCFF507QES.tar.gz
#cd ..

#mkdir -p liver2
#cd liver2
#wget https://www.encodeproject.org/files/ENCFF133HRL/@@download/ENCFF133HRL.tar.gz
#tar -xzvf ENCFF133HRL.tar.gz
#cd ..

mkdir -p liver3
#cd liver3
#wget https://www.encodeproject.org/files/ENCFF858PYR/@@download/ENCFF858PYR.tar.gz
#tar -xzvf ENCFF858PYR.tar.gz
#cd ..

#mkdir -p liver4
#cd liver4
#wget https://www.encodeproject.org/files/ENCFF448KFC/@@download/ENCFF448KFC.tar.gz
#tar -xzvf ENCFF448KFC.tar.gz
#cd ..

#mkdir -p liver5
#cd liver5
#wget https://www.encodeproject.org/files/ENCFF872JQP/@@download/ENCFF872JQP.tar.gz
#tar -xzvf ENCFF872JQP.tar.gz
#cd ..

#mkdir -p liver6
#cd liver6
#wget https://www.encodeproject.org/files/ENCFF524EJI/@@download/ENCFF524EJI.tar.gz
#tar -xzvf ENCFF524EJI.tar.gz
#cd ..

mkdir -p liver7
#cd liver7
#wget https://www.encodeproject.org/files/ENCFF728TIV/@@download/ENCFF728TIV.tar.gz
#tar -xzvf ENCFF728TIV.tar.gz
#cd ..

#mkdir -p liver8
#cd liver8
#wget https://www.encodeproject.org/files/ENCFF549XVL/@@download/ENCFF549XVL.tar.gz
#tar -xzvf ENCFF549XVL.tar.gz
#cd ..

#mkdir -p liver9
#cd liver9
#wget https://www.encodeproject.org/files/ENCFF589LJM/@@download/ENCFF589LJM.tar.gz
#tar -xzvf ENCFF589LJM.tar.gz
#cd ..

mkdir -p liver10
#cd liver10
#wget https://www.encodeproject.org/files/ENCFF429WMQ/@@download/ENCFF429WMQ.tar.gz
#tar -xzvf ENCFF429WMQ.tar.gz
#cd ..

# liver
#sample_IDs=("ENCSR594BJF-1" "ENCSR659ANG-1" "ENCSR074DOR-1" "ENCSR497ZST-1" "ENCSR825WYQ-1" 
#            "ENCSR118PUH-1" "ENCSR650BBI-1" "ENCSR552QVU-1" "ENCSR630YEA-1" "ENCSR139ATR-1")
#sample_n=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
#cell_type="liver"