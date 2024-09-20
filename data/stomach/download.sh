cd ~/scratch/scATAC_celltypes/data/stomach

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF758MUT/@@download/ENCFF758MUT.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF103XHP/@@download/ENCFF103XHP.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF206MHF/@@download/ENCFF206MHF.tar.gz #3
wget https://www.encodeproject.org/files/ENCFF423EAR/@@download/ENCFF423EAR.tar.gz #4
tar -xzvf ENCFF758MUT.tar.gz
tar -xzvf ENCFF103XHP.tar.gz
tar -xzvf ENCFF206MHF.tar.gz
tar -xzvf ENCFF423EAR.tar.gz

mkdir -p stomach1
mkdir -p stomach2
mkdir -p stomach3
mkdir -p stomach4

# stomach
#sample_IDs=("ENCSR848EID" "ENCSR549JEQ" "ENCSR266ZPZ" "ENCSR940QRM")
#sample_n=("1" "2" "3" "4")
#cell_type="stomach"