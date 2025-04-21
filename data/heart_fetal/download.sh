# download heart_fetal

mkdir -p heart_fetal1
cd heart_fetal1
wget https://www.encodeproject.org/files/ENCFF401MVS/@@download/ENCFF401MVS.tar.gz # atac
tar -xzvf ENCFF401MVS.tar.gz
mv results/ENCSR515SNH-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv results/ENCSR515SNH-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf results
rm ENCFF401MVS.tar.gz

wget https://www.encodeproject.org/files/ENCFF887HJU/@@download/ENCFF887HJU.tar.gz # rna
tar -xzvf ENCFF887HJU.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF887HJU.tar.gz
cd ..

mkdir -p heart_fetal3
cd heart_fetal3
wget https://www.encodeproject.org/files/ENCFF904ARJ/@@download/ENCFF904ARJ.tar.gz # atac
tar -xzvf ENCFF904ARJ.tar.gz
mv results/ENCSR282FAK-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv results/ENCSR282FAK-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf results
rm ENCFF904ARJ.tar.gz

wget https://www.encodeproject.org/files/ENCFF878AEW/@@download/ENCFF878AEW.tar.gz # rna
tar -xzvf ENCFF878AEW.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF878AEW.tar.gz
cd ..

mkdir -p heart_fetal7
cd heart_fetal7
wget https://www.encodeproject.org/files/ENCFF958UTA/@@download/ENCFF958UTA.tar.gz # atac
tar -xzvf ENCFF958UTA.tar.gz
mv results/ENCSR376IBI-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv results/ENCSR376IBI-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf results
rm ENCFF958UTA.tar.gz

wget https://www.encodeproject.org/files/ENCFF388YNA/@@download/ENCFF388YNA.tar.gz # rna
tar -xzvf ENCFF388YNA.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF388YNA.tar.gz
cd ..

mkdir -p heart_fetal9
cd heart_fetal9
wget https://www.encodeproject.org/files/ENCFF109FJX/@@download/ENCFF109FJX.tar.gz # atac
tar -xzvf ENCFF109FJX.tar.gz
mv results/ENCSR805DID-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv results/ENCSR805DID-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf results
rm ENCFF109FJX.tar.gz

wget https://www.encodeproject.org/files/ENCFF930KQD/@@download/ENCFF930KQD.tar.gz # rna
tar -xzvf ENCFF930KQD.tar.gz
mv GeneFull_Ex50pAS/filtered/UniqueAndMult-EM.mtx gex.mtx
mv GeneFull_Ex50pAS/filtered/barcodes.tsv rna_cells.txt
mv GeneFull_Ex50pAS/filtered/features.tsv genes.tsv
rm -rf GeneFull_Ex50pAS
rm ENCFF930KQD.tar.gz
cd ..
