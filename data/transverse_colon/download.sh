cd ~/scratch/scATAC_celltypes/data/transverse_colon

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF888LMJ/@@download/ENCFF888LMJ.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF490PBS/@@download/ENCFF490PBS.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF389FMC/@@download/ENCFF389FMC.tar.gz #3
wget https://www.encodeproject.org/files/ENCFF686ZFF/@@download/ENCFF686ZFF.tar.gz #4
wget https://www.encodeproject.org/files/ENCFF551UBJ/@@download/ENCFF551UBJ.tar.gz #5
tar -xzvf ENCFF888LMJ.tar.gz
tar -xzvf ENCFF490PBS.tar.gz
tar -xzvf ENCFF389FMC.tar.gz
tar -xzvf ENCFF686ZFF.tar.gz
tar -xzvf ENCFF551UBJ.tar.gz

mkdir -p transverse_colon1
mkdir -p transverse_colon2
mkdir -p transverse_colon3
mkdir -p transverse_colon4
mkdir -p transverse_colon5

# transverse colon
#sample_IDs=("ENCSR434SXE" "ENCSR997YNO" "ENCSR506YMX" "ENCSR007QIO" "ENCSR349XKD")
#sample_n=("1" "2" "3" "4" "5")
#cell_type="transverse_colon"