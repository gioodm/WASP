#!/usr/bin/env bash
# @Author: Giorgia Del Missier

####---- SETTING HELP MESSAGE ----####

show_help() {
cat << EOF
Usage: ./${0##*/} [-h] [-e evalue_threshold] [-b bitscore_threshold] [-t tmscore] taxid gaps_file

This script uses foldseek structural alignment results to perform gap-filling in the Genome-scale Metabolic model of interest.

    -h                      display this help and exit
    -e evalue_threshold     set the evalue threshold (default: 10e-05)
    -b bitscore_threshold   set the bitscore threshold (default: 50)
    -t tmscore_threshold    set the tmscore threshold (default: 0.5)
    taxid                   specify the GEM taxid (required)
    gaps_file               .tsv file containing the orphan reactions identified in the GEM (required)
                            --> required columns:   #rxn     rxn (extended name)            Rhea ID                                          EC number
                                             e.g:   r_0070   4-hydroxybenzoate formation    RHEA:11948, RHEA:11950, RHEA:11951, RHEA:11949   EC:3.1.2.23

Examples:
    chmod +x ${0##*/}

    ./${0##*/} 559292 559292_gaps.txt
    ./${0##*/} -e 1e-50 -b 1000 -t 0.5 559292 559292_gaps.txt
    ./${0##*/} -t 0.8 559292 559292_gaps.txt

EOF
}

export PATH=$(pwd)/foldseek/bin/:$PATH 

####---- CHECK REQUIRED COMMANDS ----####

command -v foldseek >/dev/null 2>&1 || { echo >&2 "foldseek is required but it's not installed. Aborting."; exit 1; }
command -v gsutil >/dev/null 2>&1 || { echo >&2 "gsutil is required but it's not installed. Aborting."; exit 1; }

####---- SETTING VARIABLES ----####

# Display help if -h option is provided
if [[ "$*" == *-h* ]]; then
    show_help
    exit 0
fi

# Set default values for non-positional arguments
evalue_threshold=10e-05
bitscore_threshold=50
tmscore_threshold=0.5

# Parse command-line arguments
while getopts ":e:b:t:" opt; do
    case $opt in
        e) evalue_threshold="$OPTARG";;
        b) bitscore_threshold="$OPTARG";;
        t) tmscore_threshold="$OPTARG";;
        \?) echo "Invalid option -$OPTARG" >&2; show_help; exit 1;;
    esac
done

# Remove parsed options from the positional parameters
shift $((OPTIND-1))

# Check if all required arguments are provided
if [ $# -ne 2 ]; then
    echo "Invalid option" >&2; show_help; exit 1
fi

# Assign positional arguments to variables
taxid=$1
gaps_file=$2

echo "Setting required variables:"
echo ""
echo "Selected taxid: $taxid"
echo "File containing orphan reactions: $gaps_file"
echo "Selected evalue threshold: $evalue_threshold"
echo "Selected bitscore threshold: $bitscore_threshold"
echo "Selected tmscore threshold: $tmscore_threshold"

####---- DOWNLOADING FILES AND DATABASES ----####

echo ""
echo "Downloading required foldseek databases and AlphaFold models:"
echo ""

# Define directory names
db_dir="foldseek_dbs"
prot_dir="proteomes"
results_dir="results"

# Create directories if they don't exist
mkdir -p "$db_dir" "$prot_dir" "$results_dir"

# Download foldseek databases if not already present
if [ ! -f "$db_dir/swissprot" ]; then
    foldseek databases Alphafold/Swiss-Prot "$db_dir/swissprot" tmp --remove-tmp-files 1
    rm -rf tmp
else
    echo "Foldseek databases already downloaded"
fi

# Download AlphaFold models if not already present
if [ ! -f "$prot_dir/$taxid.tar" ] && [ ! -f "$db_dir/$taxid" ]; then
    mkdir -p "$prot_dir/$taxid"
    gsutil -m cp gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-"$taxid-*"_v4.tar $prot_dir/
    for f in "$prot_dir/proteome-tax_id-$taxid"-*_v4.tar; do
        tar -xf "$f" -C "$prot_dir/$taxid"
        rm "$f"
    done
    for f in "$prot_dir/$taxid"/*.json.gz; do
        rm "$f" 
    done
    tar -cf "$prot_dir/$taxid.tar" -C "$prot_dir/" "$taxid"
    rm -r "$prot_dir/$taxid"

    foldseek createdb "$prot_dir/$taxid.tar" "$db_dir/$taxid"
elif [ -f "$prot_dir/$taxid.tar" ] && [ ! -f "$db_dir/$taxid" ]; then
    foldseek createdb "$prot_dir/$taxid.tar" "$db_dir/$taxid"
else
    echo "AlphaFold models of selected organism already downloaded"
fi

####---- UNIPROT MAPPING ----####

echo ""
echo "Mapping orphan reaction to UniProt IDs"
echo ""

taxid_dir="$results_dir/$taxid"
mkdir -p "$taxid_dir"

python3 rxn2uniprot.py --input "$gaps_file" --output "$taxid_dir/${taxid}_rxn2up.txt" --output_ids "$taxid_dir/${taxid}_upIDs.txt"

echo "done"

####---- FOLDSEEK STRUCTURAL ALIGNMENT ----####

echo ""
echo "Performing foldseek alignment search:"
echo ""

foldseek prefixid "$db_dir/swissprot"_h "$db_dir/swissprot.lookup" --tsv --threads 1
awk 'NR == FNR {f[$1] = $1; next} $2 in f {print $1}' "$taxid_dir/${taxid}_upIDs.txt" "$db_dir/swissprot.lookup" > "$db_dir/subset$taxid.tsv"

foldseek createsubdb "$db_dir/subset$taxid.tsv" "$db_dir/swissprot" "$db_dir/subdb$taxid"
foldseek createsubdb "$db_dir/subset$taxid.tsv" "$db_dir/swissprot_ss" "$db_dir/subdb$taxid_ss"
foldseek createsubdb "$db_dir/subset$taxid.tsv" "$db_dir/swissprot_ca" "$db_dir/subdb$taxid_ca"
rm "$db_dir/subset$taxid.tsv"

foldseek search "$db_dir/subdb$taxid" "$db_dir/$taxid" "$taxid_dir/$taxid" tmp -a 1 --threads 64
foldseek convertalis --format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits" \
                     "$db_dir/subdb$taxid" "$db_dir/$taxid" "$taxid_dir/$taxid" "$taxid_dir/$taxid.m8"

foldseek search "$db_dir/subdb$taxid" "$db_dir/subdb$taxid" "$taxid_dir/${taxid}_db_allvsall" tmp -a 1 --threads 64
foldseek convertalis --format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits" \
                     "$db_dir/subdb$taxid" "$db_dir/subdb$taxid" "$taxid_dir/${taxid}_db_allvsall" "$taxid_dir/${taxid}_db_allvsall.m8"
    
rm "$db_dir/subdb$taxid"*
rm "$db_dir/$taxid"*
find "$taxid_dir" -type f ! -name "*.txt" ! -name "*.m8" -delete

####---- GAP FILLING ----####

echo ""
echo "Identifying best hits for each orphan reaction and performing gap-filling"
echo ""

python3 gap_filling.py --input "$taxid_dir/$taxid.m8" --input_db "$taxid_dir/${taxid}_db_allvsall.m8" --input_rxn2up "$taxid_dir/${taxid}_rxn2up.txt" \
                       --output "$taxid_dir/${taxid}_hits.txt" --eval_thr $evalue_threshold --bits_thr $bitscore_threshold --tm_thr $tmscore_threshold

echo "done"
echo ""
echo "GEM gap-filling complete!"
