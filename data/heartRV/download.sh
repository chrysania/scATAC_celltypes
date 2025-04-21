# download heartRV

mkdir -p heartRV1
cd heartRV1
wget https://www.encodeproject.org/files/ENCFF223DBM/@@download/ENCFF223DBM.tar.gz # atac
tar -xzvf ENCFF223DBM.tar.gz
mv encode_scatac_dcc_2/results/ENCSR588PEE-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR588PEE-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF223DBM.tar.gz

wget https://www.encodeproject.org/files/ENCFF757DCI/@@download/ENCFF757DCI.tar.gz # rna
tar -xzvf ENCFF757DCI.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF757DCI.tar.gz
cd .. 

mkdir -p heartRV2
cd heartRV2
wget https://www.encodeproject.org/files/ENCFF145ABZ/@@download/ENCFF145ABZ.tar.gz # atac
tar -xzvf ENCFF145ABZ.tar.gz
mv encode_scatac_dcc_2/results/ENCSR681OLJ-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR681OLJ-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF145ABZ.tar.gz

wget https://www.encodeproject.org/files/ENCFF289QCB/@@download/ENCFF289QCB.tar.gz # rna
tar -xzvf ENCFF289QCB.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF289QCB.tar.gz
cd .. 

mkdir -p heartRV5
cd heartRV5
wget https://www.encodeproject.org/files/ENCFF358RDK/@@download/ENCFF358RDK.tar.gz # atac
tar -xzvf ENCFF358RDK.tar.gz
mv encode_scatac_dcc_2/results/ENCSR814OLA-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR814OLA-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF358RDK.tar.gz
cd ..

mkdir -p heartRV9
cd heartRV9
wget https://www.encodeproject.org/files/ENCFF655IKC/@@download/ENCFF655IKC.tar.gz # atac
tar -xzvf ENCFF655IKC.tar.gz
mv encode_scatac_dcc_2/results/ENCSR517QNQ-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR517QNQ-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF655IKC.tar.gz
cd ..