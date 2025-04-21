# download liver

mkdir -p liver3
cd liver3
wget https://www.encodeproject.org/files/ENCFF589YDE/@@download/ENCFF589YDE.tar.gz # atac
tar -xzvf ENCFF589YDE.tar.gz
mv encode_scatac_dcc_2/results/ENCSR074DOR-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR074DOR-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF589YDE.tar.gz

wget https://www.encodeproject.org/files/ENCFF858PYR/@@download/ENCFF858PYR.tar.gz # rna
tar -xzvf ENCFF858PYR.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF858PYR.tar.gz
cd ..

mkdir -p liver7
cd liver7
wget https://www.encodeproject.org/files/ENCFF883DOK/@@download/ENCFF883DOK.tar.gz # atac
tar -xzvf ENCFF883DOK.tar.gz
mv encode_scatac_dcc_2/results/ENCSR650BBI-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR650BBI-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF883DOK.tar.gz

wget https://www.encodeproject.org/files/ENCFF728TIV/@@download/ENCFF728TIV.tar.gz # rna
tar -xzvf ENCFF728TIV.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF728TIV.tar.gz
cd ..