# download adrenal

mkdir -p adrenal1
cd adrenal1
wget https://www.encodeproject.org/files/ENCFF318DQR/@@download/ENCFF318DQR.tar.gz # atac
tar -xzvf ENCFF318DQR.tar.gz
mv encode_scatac_dcc_2/results/ENCSR420EWQ-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR420EWQ-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm -rf ENCFF318DQR.tar.gz

wget https://www.encodeproject.org/files/ENCFF451NNZ/@@download/ENCFF451NNZ.tar.gz # rna
tar -xzvf ENCFF451NNZ.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF451NNZ.tar.gz
cd ..

mkdir -p adrenal2
cd adrenal2
wget https://www.encodeproject.org/files/ENCFF370JZO/@@download/ENCFF370JZO.tar.gz # atac
tar -xzvf ENCFF370JZO.tar.gz
mv encode_scatac_dcc_2/results/ENCSR693GAD-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR693GAD-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF370JZO.tar.gz

wget https://www.encodeproject.org/files/ENCFF846BGR/@@download/ENCFF846BGR.tar.gz # rna
tar -xzvf ENCFF846BGR.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm ENCFF846BGR.tar.gz
cd ..