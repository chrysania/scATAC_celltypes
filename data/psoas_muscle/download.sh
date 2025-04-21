# download psoas_muscle

mkdir -p psoas_muscle1
cd psoas_muscle1
wget https://www.encodeproject.org/files/ENCFF948PTQ/@@download/ENCFF948PTQ.tar.gz # atac
tar -xzvf ENCFF948PTQ.tar.gz
mv encode_scatac_dcc_2/results/ENCSR000XQD-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR000XQD-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF948PTQ.tar.gz

wget https://www.encodeproject.org/files/ENCFF644KGC/@@download/ENCFF644KGC.tar.gz # rna
tar -xzvf ENCFF644KGC.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF644KGC.tar.gz
cd ..

mkdir -p psoas_muscle2
cd psoas_muscle2
wget https://www.encodeproject.org/files/ENCFF769KLR/@@download/ENCFF769KLR.tar.gz # atac
tar -xzvf ENCFF769KLR.tar.gz
mv encode_scatac_dcc_2/results/ENCSR869LEV-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR869LEV-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF769KLR.tar.gz

wget https://www.encodeproject.org/files/ENCFF810VOM/@@download/ENCFF810VOM.tar.gz # rna
tar -xzvf ENCFF810VOM.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF810VOM.tar.gz
cd ..

mkdir -p psoas_muscle3
cd psoas_muscle3
wget https://www.encodeproject.org/files/ENCFF832YIV/@@download/ENCFF832YIV.tar.gz # atac
tar -xzvf ENCFF832YIV.tar.gz
mv encode_scatac_dcc_2/results/ENCSR916EDP-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR916EDP-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF832YIV.tar.gz

wget https://www.encodeproject.org/files/ENCFF831OBZ/@@download/ENCFF831OBZ.tar.gz # rna
tar -xzvf ENCFF831OBZ.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF831OBZ.tar.gz
cd ..