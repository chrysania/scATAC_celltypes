cd ~/scratch/scATAC_celltypes/data/peyer

# ATAC fragments
wget https://www.encodeproject.org/files/ENCFF865UKG/@@download/ENCFF865UKG.tar.gz #1
wget https://www.encodeproject.org/files/ENCFF286BTT/@@download/ENCFF286BTT.tar.gz #2
wget https://www.encodeproject.org/files/ENCFF731NML/@@download/ENCFF731NML.tar.gz #3
tar -xzvf ENCFF865UKG.tar.gz
tar -xzvf ENCFF286BTT.tar.gz
tar -xzvf ENCFF731NML.tar.gz

mkdir -p peyer1
mkdir -p peyer2
mkdir -p peyer3

# peyer's patch
#sample_IDs=("ENCSR052DKH" "ENCSR652WJE" "ENCSR101JHK")
#sample_n=("1" "2" "3")
#cell_type="peyer"