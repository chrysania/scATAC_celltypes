# move to code dir
cd ~/scratch/scATAC_celltypes/code

# function to call peaks and build matrix
callpeaks() {

    cell_type=$1
    sample_id=$2
    sample_n=$3
    
    cd ../data/$cell_type/
    pwd

    echo "dataset ID = $sample_id"
    
    n=${sample_n[i]}
    echo "dataset number = $n"
    
    frag_path="encode_scatac_dcc_2/results/${sample_id}/fragments/fragments.tsv.gz"
    echo "fragment path = $frag_path"

    feature_path="$cell_type${n}/$cell_type${n}_peakcalling.txt"
    echo $feature_path

    echo "running fragtk count on dataset $sample_id"
    ~/scratch/fragtk/target/release/fragtk count \
        -f $frag_path \
        -o $cell_type${n}/barcode_counts.tsv > /dev/null

    echo "getting top 10,000 barcodes"
    sort -k2,2nr $cell_type${n}/barcode_counts.tsv | head -n 10000 > $cell_type${n}/10000_barcode_counts.tsv
    cut -f1 $cell_type${n}/10000_barcode_counts.tsv > $cell_type${n}/10000_barcodes_only.tsv

    echo "calling peaks with macs2"
    macs2 callpeak \
    -f BED --nomodel --shift -100 --extsize 200 --name $cell_type${n} \
    -t $frag_path \
	--outdir $cell_type${n}/
    
    cut -f1-3 $cell_type${n}/$cell_type${n}_peaks.narrowPeak > $cell_type${n}/$cell_type${n}_peakcalling.txt

    # fragtk matrix to generate peak by cell matrix
    echo "running fragtk matrix"
    ~/scratch/fragtk/target/release/fragtk matrix \
        -f $frag_path \
        -b $feature_path \
        -c $cell_type${n}/10000_barcodes_only.tsv \
        -o $cell_type${n}/peak_matrix

    cd ../../code/
    # make UMAP
    echo "making UMAP plot"
    Rscript celltype_umap.r ${n} $cell_type
}

# esophagus
sample_IDs=("ENCSR453TVZ-1" "ENCSR757EGB-1" "ENCSR164GSH-1")
sample_n=("1" "2" "3")
cell_type="esophagus"

for i in $(seq 0 2); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# heart_fetal 
sample_IDs=("ENCSR515SNH-1" "ENCSR715JSZ-1" "ENCSR282FAK-1" "ENCSR890TGR-1" "ENCSR024TGD-1" "ENCSR004IAY-1" "ENCSR376IBI-1" "ENCSR306FRQ-1" "ENCSR805DID-1")
sample_n=("1" "2" "3" "4" "5" "6" "7" "8" "9")
cell_type="heart_fetal"

for i in $(seq 0 8); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# heartLV (left ventricle)
sample_IDs=("ENCSR321AHR-1" "ENCSR769WLL-1" "ENCSR701JAT-1" "ENCSR270CUT-1" "ENCSR960IDI-1" "ENCSR506ROZ-1" "ENCSR862IWS-1" "ENCSR088ZOL-1" "ENCSR627IOJ-1" "ENCSR020FAW-1")
sample_n=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
cell_type="heartLV"

for i in $(seq 0 9); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# heartRV (right ventricle)
sample_IDs=("ENCSR588PEE-1" "ENCSR681OLJ-1" "ENCSR169BCG-1" "ENCSR604PDO-1" "ENCSR814OLA-1" "ENCSR520ZUD-1" "ENCSR579TPC-1" "ENCSR615TSN-1" "ENCSR517QNQ-1" "ENCSR454YDZ-1")
sample_n=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
cell_type="heartRV"

for i in $(seq 0 9); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# left colon
sample_IDs=("ENCSR830FPR-1" "ENCSR916RYB-1" "ENCSR904WIW-1")
sample_n=("1" "2" "3")
cell_type="left_colon"

for i in $(seq 0 2); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# lower leg skin
sample_IDs=("ENCSR733SZL-1" "ENCSR397ODX-1" "ENCSR474TGL-1" "ENCSR513HZN-1")
sample_n=("1" "2" "3" "4")
cell_type="leg_skin"

for i in $(seq 0 3); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# gastrocnemius medialis (muscle2)
sample_IDs=("ENCSR696YOC-1" "ENCSR819EGE-1" "ENCSR244GZL-1" "ENCSR023FME-1" "ENCSR139TIQ-1")
sample_n=("1" "2" "3" "4" "5")
cell_type="leg_muscle"

for i in $(seq 0 4); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# omental fat pad
sample_IDs=("ENCSR181XXQ-1" "ENCSR492GGN-1" "ENCSR644SCP-1" "ENCSR274HQD-1")
sample_n=("1" "2" "3" "4")
cell_type="omental_fat"

for i in $(seq 0 3); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# peyer's patch
sample_IDs=("ENCSR052DKH-1" "ENCSR652WJE-1" "ENCSR101JHK-1")
sample_n=("1" "2" "3")
cell_type="peyer"

for i in $(seq 0 2); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# stomach
sample_IDs=("ENCSR848EID-1" "ENCSR549JEQ-1" "ENCSR266ZPZ-1" "ENCSR940QRM-1")
sample_n=("1" "2" "3" "4")
cell_type="stomach"

for i in $(seq 0 3); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# thyroid
sample_IDs=("ENCSR796RXX-1" "ENCSR909OXO-1" "ENCSR817VFO-1")
sample_n=("1" "2" "3")
cell_type="thyroid"

for i in $(seq 0 2); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# tibial nerve
sample_IDs=("ENCSR205TUH-1" "ENCSR453BVR-1" "ENCSR726QTF-1")
sample_n=("1" "2" "3")
cell_type="tibial_nerve"

for i in $(seq 0 2); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# transverse colon
sample_IDs=("ENCSR434SXE-1" "ENCSR997YNO-1" "ENCSR506YMX-1" "ENCSR007QIO-1" "ENCSR349XKD-1")
sample_n=("1" "2" "3" "4" "5")
cell_type="transverse_colon"

for i in $(seq 0 4); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# uterus
sample_IDs=("ENCSR828HVB-1" "ENCSR455CVZ-1" "ENCSR028RFK-1")
sample_n=("1" "2" "3")
cell_type="uterus"

for i in $(seq 0 2); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

############################# from metapeak-analysis
# adrenal gland
sample_IDs=("ENCSR420EWQ-1" "ENCSR693GAD-1" "ENCSR194KHA-1")
sample_n=("1" "2" "3")
cell_type="adrenal"

for i in $(seq 0 2); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# lung
sample_IDs=("ENCSR816NWE-1" "ENCSR391BWM-1" "ENCSR824OCY-1")
sample_n=("1" "2" "3")
cell_type="lung"

for i in $(seq 0 2); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# psoas muscle (muscle1)
sample_IDs=("ENCSR000XQD-1" "ENCSR869LEV-1" "ENCSR916EDP-1" "ENCSR332XEW-1")
sample_n=("1" "2" "3" "4")
cell_type="psoas_muscle"

for i in $(seq 0 3); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# ovary
#sample_IDs=("ENCSR322LEY-1" "ENCSR533CSY-1" "ENCSR305VES-1" "ENCSR105VKN-1" "ENCSR197GLE-1")
#sample_n=("1" "2" "3" "4" "5")
#cell_type="ovary"

#for i in $(seq 0 4); do
#    sample_id=${sample_IDs[i]}
#    callpeaks $cell_type $sample_id    
#done

# pancreas
sample_IDs=("ENCSR868CRK-1" "ENCSR690NZI-1" "ENCSR229VVY-1" "ENCSR496XXB-1")
sample_n=("1" "2" "3" "4")
cell_type="pancreas"

for i in $(seq 0 3); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done

# liver
sample_IDs=("ENCSR594BJF-1" "ENCSR659ANG-1" "ENCSR074DOR-1" "ENCSR497ZST-1" "ENCSR825WYQ-1" 
            "ENCSR118PUH-1" "ENCSR650BBI-1" "ENCSR552QVU-1" "ENCSR630YEA-1" "ENCSR139ATR-1")
sample_n=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
cell_type="liver"

for i in $(seq 0 9); do
    sample_id=${sample_IDs[i]}
    callpeaks $cell_type $sample_id    
done
