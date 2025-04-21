# download esophagus

mkdir -p esophagus2
cd esophagus2
wget https://www.encodeproject.org/files/ENCFF623PSE/@@download/ENCFF623PSE.tar.gz
tar -xzvf ENCFF623PSE.tar.gz
mv encode_scatac_dcc_2/results/ENCSR757EGB-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR757EGB-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF623PSE.tar.gz
cd ..

mkdir -p esophagus3
cd esophagus3
wget https://www.encodeproject.org/files/ENCFF815ESF/@@download/ENCFF815ESF.tar.gz
tar -xzvf ENCFF815ESF.tar.gz
mv encode_scatac_dcc_2/results/ENCSR164GSH-1/fragments/fragments.tsv.gz fragments.tsv.gz
mv encode_scatac_dcc_2/results/ENCSR164GSH-1/fragments/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
rm -rf encode_scatac_dcc_2
rm ENCFF815ESF.tar.gz
cd ..