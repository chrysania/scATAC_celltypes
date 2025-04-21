# download left colon

mkdir -p left_colon1
cd left_colon1
wget https://www.encodeproject.org/files/ENCFF040VYY/@@download/ENCFF040VYY.tar.gz # atac
tar -xzvf ENCFF040VYY.tar.gz
mv encode_scatac_dcc_2/results/ENCSR830FPR-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR830FPR-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF040VYY.tar.gz

wget https://www.encodeproject.org/files/ENCFF949IOU/@@download/ENCFF949IOU.tar.gz # rna
tar -xzvf ENCFF949IOU.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF949IOU.tar.gz
cd ..

mkdir -p left_colon2
wget https://www.encodeproject.org/files/ENCFF826JKT/@@download/ENCFF826JKT.tar.gz # atac
tar -xzvf ENCFF826JKT.tar.gz
mv encode_scatac_dcc_2/results/ENCSR916RYB-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR916RYB-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF826JKT.tar.gz
cd ..