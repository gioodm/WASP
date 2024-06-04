#!/usr/bin/env bash
# @author Giorgia Del Missier


####---- SETTING HELP MESSAGE ----####

show_help() {
cat << EOF
Usage: ./${0##*/} [-h] [-e evalue_threshold] [-b bitscore_threshold] [-t tmscore] taxid gaps_file

This script uses foldseek structural alignment results to perform gap-filling in a GEnome-scale Metabolic model of interest.
It also includes optional evalue/bitscore/tmscore thresholds.

    -h                      display this help and exit
    -e evalue_threshold     set the evalue threshold (default: 10e-05)
    -b bitscore_threshold   set the bitscore threshold (default: 50)
    -t tmscore_threshold    set the tmscore threshold (default: 0.5)
    taxid                   specify the GEM taxid (required)
    gaps_file               .txt file containing the orphan reactions identified in the GEM (required)
                            --> required columns:   #rxn     rxn (extended name)            Rhea ID                                          EC number
                                             e.g:   r_0070   4-hydroxybenzoate formation    RHEA:11948, RHEA:11950, RHEA:11951, RHEA:11949   EC:3.1.2.23

Examples:
    chmod +x ${0##*/}

    ./${0##*/} 559292 559292_gaps.txt
    ./${0##*/} -e 1e-50 -b 1000 -t 0.5 559292 559292_gaps.txt
    ./${0##*/} -t 0.8 559292 559292_gaps.txt

EOF
}


####---- SETTING VARIABLES ----####

# check if -h option is provided
if [[ "$*" == *-h* ]]; then
    show_help
    exit 0
fi

# set default values for non-positional arguments (evalue threshold, bitscore threshold and tmscore threshold)
e=10e-05
b=50
t=0.5

# parse command-line arguments
while getopts ":b:e:t:" opt; do
    case $opt in
        e) e="$OPTARG";;

        b) b="$OPTARG";;

        t) t="$OPTARG";;

        \?) echo "Invalid option -$OPTARG" >&2; show_help; exit 1;;
    esac
done

# remove the parsed options from the positional parameters
shift $((OPTIND-1))

# check if all required arguments are provided
if [ $# -ne 2 ]; then
    echo "Invalid option" >&2; show_help; exit 1
fi

# assign positional arguments to variables
taxid=$1
gaps_file=$2

echo "Setting required variables:"
echo ""
echo "Selected taxid is: $taxid"
echo "File containig orphan reactions set to: $gaps_file"
echo "Selected evalue threshold is: $e"
echo "Selected bitscore threshold is: $b"
echo "Selected tmscore threshold is: $t"


####---- DOWNLOADING FILES AND DATABASES ----####

echo ""
echo "Downloading required foldseek databases and AlphaFold models:"
echo ""

export PATH=$(pwd)/foldseek/bin/:$PATH 

# define directory names
db_dir="foldseek_dbs"
prot_dir="proteomes"
results_dir="results"

# check if directories exist, create if necessary
if [ ! -d "$db_dir" ]; then
    mkdir "$db_dir"
fi

if [ ! -d "$prot_dir" ]; then
    mkdir "$prot_dir"
fi

if [ ! -d "$results_dir" ]; then
    mkdir "$results_dir"
fi

# check if files exist, download if necessary
if [ ! -f "$db_dir/swissprot" ]; then
    foldseek databases Alphafold/Swiss-Prot "$db_dir/swissprot" tmp
    rm -r tmp
else
    echo "foldseek databases already downloaded"
fi

if [ ! -f "$prot_dir/$taxid.tar" ] && [ ! -f "$db_dir/$taxid" ]; then
    mkdir "$prot_dir/$taxid"
    gsutil -m cp gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-"$taxid"-*_v4.tar $prot_dir/
    for f in "$prot_dir/proteome-tax_id-$taxid"-*_v4.tar; do
        tar -xf "$f" -C "$prot_dir/$taxid"
        rm "$f"
    done
    for f in "$prot_dir/$taxid"/*.json.gz; do
        rm "$f" 
    done
    tar -cf "$prot_dir/$taxid.tar" -C "$prot_dir/" "$taxid"
    rm -r "$prot_dir/$taxid"

    foldseek createdb "$prot_dir/$taxid.tar" "$db_dir/$taxid" --tar-include ".*cif"

elif [ -f "$prot_dir/$taxid.tar" ] && [ ! -f "$db_dir/$taxid" ]; then

    foldseek createdb "$prot_dir/$taxid.tar" "$db_dir/$taxid" --tar-include ".*cif"

else
    echo "AlphaFold models of selected organism already downloaded"
fi


####---- UNIPROT MAPPING ----####

echo ""
echo "Mapping orphan reaction to UniProt IDS"
echo ""

taxid_dir="$results_dir/$taxid"

if [ ! -d "$taxid_dir" ]; then
    mkdir "$taxid_dir"
fi

python3 rxn2uniprot.py --input $gaps_file --output "$taxid_dir/$taxid"_rxn2up.txt --output_ids "$taxid_dir/$taxid"_upIDs.txt

echo "done"


####---- FOLDSEEK STRUCTURAL ALIGNMENT ----####

echo ""
echo "Performing foldseek alignment search:"
echo ""

foldseek prefixid "$db_dir/swissprot"_h "$db_dir/swissprot".lookup --tsv --threads 1
awk 'NR == FNR {f[$1] = $1; next} $2 in f {print $1}' "$taxid_dir/$taxid"_upIDs.txt "$db_dir/swissprot".lookup > "$db_dir/subset$taxid".tsv
foldseek createsubdb "$db_dir/subset$taxid".tsv "$db_dir/swissprot" "$db_dir/subdb$taxid"
foldseek createsubdb "$db_dir/subset$taxid".tsv "$db_dir/swissprot"_ss "$db_dir/subdb$taxid"_ss
foldseek createsubdb "$db_dir/subset$taxid".tsv "$db_dir/swissprot"_ca "$db_dir/subdb$taxid"_ca
rm "$db_dir/subset$taxid".tsv

foldseek search "$db_dir/subdb$taxid" "$db_dir/$taxid" "$taxid_dir/$taxid" tmp -a 1 --threads 64
foldseek convertalis --format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits" \
                     "$db_dir/subdb$taxid" "$db_dir/$taxid" "$taxid_dir/$taxid" "$taxid_dir/$taxid".m8

foldseek search "$db_dir/subdb$taxid" "$db_dir/subdb$taxid" "$taxid_dir/$taxid"_db_allvsall tmp -a 1 --threads 64
foldseek convertalis --format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits" \
                     "$db_dir/subdb$taxid" "$db_dir/subdb$taxid" "$taxid_dir/$taxid"_db_allvsall "$taxid_dir/$taxid"_db_allvsall.m8
    
rm "$db_dir/subdb$taxid"*
rm "$db_dir/$taxid"*
find "$taxid_dir" -type f ! -name "*.txt" ! -name "*.m8" -delete


####---- GAP FILLING ----####

echo ""
echo "Identifying best hits for each orphan reaction to perform gap-filling"
echo ""

python3 gap_filling.py --input "$taxid_dir/$taxid".m8 --input_db "$taxid_dir/$taxid"_db_allvsall.m8 --input_rxn2up "$taxid_dir/$taxid"_rxn2up.txt  \
                       --output "$taxid_dir/$taxid"_hits.txt --eval_thr $e --bits_thr $b --tm_thr $t

echo "done"


####---- ANNOTATION ----####

echo ""
echo "Annotating best hits"
echo ""

python3 get_annotation_GEM.py --input "$taxid_dir/$taxid"_hits.txt --output "$taxid_dir/$taxid"_hits_annotated.txt

echo "done"

