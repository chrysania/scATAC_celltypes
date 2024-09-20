cd ~/scratch/scATAC_celltypes/data/heart_fetal

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF401MVS/@@download/ENCFF401MVS.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF946HNM/@@download/ENCFF946HNM.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF904ARJ/@@download/ENCFF904ARJ.tar.gz #3
wget https://www.encodeproject.org/files/ENCFF658TEH/@@download/ENCFF658TEH.tar.gz #4
wget https://www.encodeproject.org/files/ENCFF908JHS/@@download/ENCFF908JHS.tar.gz #5
wget https://www.encodeproject.org/files/ENCFF851VTB/@@download/ENCFF851VTB.tar.gz #6
wget https://www.encodeproject.org/files/ENCFF958UTA/@@download/ENCFF958UTA.tar.gz #7
wget https://www.encodeproject.org/files/ENCFF575DGZ/@@download/ENCFF575DGZ.tar.gz #8
wget https://www.encodeproject.org/files/ENCFF109FJX/@@download/ENCFF109FJX.tar.gz #9
tar -xzvf ENCFF401MVS.tar.gz
tar -xzvf ENCFF946HNM.tar.gz
tar -xzvf ENCFF904ARJ.tar.gz
tar -xzvf ENCFF658TEH.tar.gz
tar -xzvf ENCFF908JHS.tar.gz
tar -xzvf ENCFF851VTB.tar.gz
tar -xzvf ENCFF958UTA.tar.gz
tar -xzvf ENCFF575DGZ.tar.gz
tar -xzvf ENCFF109FJX.tar.gz

# RNA matrix
mkdir -p heart_fetal1
#https://www.encodeproject.org/files/ENCFF887HJU/@@download/ENCFF887HJU.tar.gz

mkdir -p heart_fetal2
#https://www.encodeproject.org/files/ENCFF243UTL/@@download/ENCFF243UTL.tar.gz

mkdir -p heart_fetal3
#https://www.encodeproject.org/files/ENCFF878AEW/@@download/ENCFF878AEW.tar.gz

mkdir -p heart_fetal4
#https://www.encodeproject.org/files/ENCFF654KZX/@@download/ENCFF654KZX.tar.gz

mkdir -p heart_fetal5
#https://www.encodeproject.org/files/ENCFF260ZOM/@@download/ENCFF260ZOM.tar.gz

mkdir -p heart_fetal6
#https://www.encodeproject.org/files/ENCFF807XXG/@@download/ENCFF807XXG.tar.gz

mkdir -p heart_fetal7
#https://www.encodeproject.org/files/ENCFF388YNA/@@download/ENCFF388YNA.tar.gz

mkdir -p heart_fetal8
#https://www.encodeproject.org/files/ENCFF795QHU/@@download/ENCFF795QHU.tar.gz

mkdir -p heart_fetal9
#https://www.encodeproject.org/files/ENCFF930KQD/@@download/ENCFF930KQD.tar.gz

# heart_fetal 
#sample_IDs=("ENCSR515SNH" "ENCSR715JSZ" "ENCSR282FAK" "ENCSR890TGR" "ENCSR024TGD" "ENCSR004IAY" "ENCSR376IBI" "ENCSR306FRQ" "ENCSR805DID")
#sample_n=("1" "2" "3" "4" "5" "6" "7" "8" "9")
#cell_type="heart_fetal"