cd ~/scratch/scATAC_celltypes/data/muscle1

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF948PTQ/@@download/ENCFF948PTQ.tar.gz
wget https://www.encodeproject.org/files/ENCFF769KLR/@@download/ENCFF769KLR.tar.gz
wget https://www.encodeproject.org/files/ENCFF832YIV/@@download/ENCFF832YIV.tar.gz
wget https://www.encodeproject.org/files/ENCFF301RZM/@@download/ENCFF301RZM.tar.gz
tar -xzvf ENCFF948PTQ.tar.gz
tar -xzvf ENCFF769KLR.tar.gz
tar -xzvf ENCFF832YIV.tar.gz
tar -xzvf ENCFF301RZM.tar.gz

# RNA matrix
mkdir -p muscle1
cd muscle1
wget https://www.encodeproject.org/files/ENCFF644KGC/@@download/ENCFF644KGC.tar.gz
tar -xzvf ENCFF644KGC.tar.gz
cd ..

mkdir -p muscle2
cd muscle2
wget https://www.encodeproject.org/files/ENCFF810VOM/@@download/ENCFF810VOM.tar.gz
tar -xzvf ENCFF810VOM.tar.gz
cd ..

mkdir -p muscle3
cd muscle3
wget https://www.encodeproject.org/files/ENCFF831OBZ/@@download/ENCFF831OBZ.tar.gz
tar -xzvf ENCFF831OBZ.tar.gz
cd ..

mkdir -p muscle4
cd muscle4
wget https://www.encodeproject.org/files/ENCFF502DDT/@@download/ENCFF502DDT.tar.gz
tar -xzvf ENCFF502DDT.tar.gz
cd ..

# psoas muscle
#sample_IDs=("ENCSR000XQD-1" "ENCSR869LEV-1" "ENCSR916EDP-1" "ENCSR332XEW-1")
#sample_n=("1" "2" "3" "4")
#cell_type="muscle1"