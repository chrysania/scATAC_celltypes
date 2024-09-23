cd ~/scratch/scATAC_celltypes/data/heartLV

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF920ZTM/@@download/ENCFF920ZTM.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF478UNG/@@download/ENCFF478UNG.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF673ECB/@@download/ENCFF673ECB.tar.gz #3
wget https://www.encodeproject.org/files/ENCFF730YBD/@@download/ENCFF730YBD.tar.gz #4
wget https://www.encodeproject.org/files/ENCFF157PSM/@@download/ENCFF157PSM.tar.gz #5
wget https://www.encodeproject.org/files/ENCFF888EOS/@@download/ENCFF888EOS.tar.gz #6
wget https://www.encodeproject.org/files/ENCFF962ZWK/@@download/ENCFF962ZWK.tar.gz #7
wget https://www.encodeproject.org/files/ENCFF296HET/@@download/ENCFF296HET.tar.gz #8
wget https://www.encodeproject.org/files/ENCFF478CZM/@@download/ENCFF478CZM.tar.gz #9
wget https://www.encodeproject.org/files/ENCFF654FKD/@@download/ENCFF654FKD.tar.gz #10
tar -xzvf ENCFF920ZTM.tar.gz
tar -xzvf ENCFF478UNG.tar.gz
tar -xzvf ENCFF673ECB.tar.gz
tar -xzvf ENCFF730YBD.tar.gz
tar -xzvf ENCFF157PSM.tar.gz
tar -xzvf ENCFF888EOS.tar.gz
tar -xzvf ENCFF962ZWK.tar.gz
tar -xzvf ENCFF296HET.tar.gz
tar -xzvf ENCFF478CZM.tar.gz
tar -xzvf ENCFF654FKD.tar.gz

# RNA matrix
mkdir -p heartLV1
#wget https://www.encodeproject.org/files/ENCFF552HTP/@@download/ENCFF552HTP.tar.gz

mkdir -p heartLV2
#https://www.encodeproject.org/files/ENCFF538RCL/@@download/ENCFF538RCL.tar.gz

mkdir -p heartLV3
#https://www.encodeproject.org/files/ENCFF625DIE/@@download/ENCFF625DIE.tar.gz

mkdir -p heartLV4
#https://www.encodeproject.org/files/ENCFF046HGL/@@download/ENCFF046HGL.tar.gz

mkdir -p heartLV5
#https://www.encodeproject.org/files/ENCFF379FXA/@@download/ENCFF379FXA.tar.gz

mkdir -p heartLV6
#https://www.encodeproject.org/files/ENCFF004QRO/@@download/ENCFF004QRO.tar.gz

mkdir -p heartLV7
#https://www.encodeproject.org/files/ENCFF323OMW/@@download/ENCFF323OMW.tar.gz

mkdir -p heartLV8
#https://www.encodeproject.org/files/ENCFF484HHL/@@download/ENCFF484HHL.tar.gz

mkdir -p heartLV9
#https://www.encodeproject.org/files/ENCFF431AQJ/@@download/ENCFF431AQJ.tar.gz

mkdir -p heartLV10
#https://www.encodeproject.org/files/ENCFF802HTG/@@download/ENCFF802HTG.tar.gz

# heartLV (left ventricle)
#sample_IDs=("ENCSR321AHR-1" "ENCSR769WLL-1" "ENCSR701JAT-1" "ENCSR270CUT-1" "ENCSR960IDI-1" "ENCSR506ROZ-1" "ENCSR862IWS-1" "ENCSR088ZOL-1" "ENCSR627IOJ-1" "ENCSR020FAW-1")
#sample_n=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
#cell_type="heartLV"