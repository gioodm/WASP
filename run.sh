#!/usr/bin/env bash
# @author Giorgia Del Missier

####---- HELP MESSAGE FUNCTION ----####
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-e evalue_threshold] [-b bitscore_threshold] [-n max_neighbours] [-s step] [-i iterations] taxid

WASP (Whole-proteome Annotation through Structural homology Pipeline) performs a "structural BLAST" using AlphaFold models to better annotate the target taxid proteome.
Parameters:

    -h                        display this help and exit
    -e evalue_threshold       set the evalue threshold (default: 10e-10)
    -b bitscore threshold     set the bitscore threshold (default: 50)
    -n max_neighbours         set the max number of neighbours (default: 10)
    -s step                   set step to add to max neighbours (n) in additional iterations (default: 10)
    -i iterations             set number of iterations to perform (default: 3)
    taxid                     NCBI taxonomy identifier to be analysed (required)

Examples:
    chmod +x ${0##*/}

    ${0##*/} 559292
    ${0##*/} -e 1e-50 -b 200 -n 5 -i 5 559292
    ${0##*/} -s 5 559292

EOF
}

export PATH=$(pwd)/foldseek/bin/:$PATH 

####---- CHECK REQUIRED COMMANDS ----####

command -v foldseek >/dev/null 2>&1 || { echo >&2 "foldseek is required but it's not installed. Aborting."; exit 1; }
command -v gsutil >/dev/null 2>&1 || { echo >&2 "gsutil is required but it's not installed. Aborting."; exit 1; }

####---- SETTING VARIABLES ----####

# Check if -h option is provided
if [[ "$*" == *-h* ]]; then
    show_help
    exit 0
fi

# Set default values for non-positional arguments
e=10e-10
b=50
n=10
s=10
i=3

# Parse command-line arguments
while getopts ":e:b:n:s:i:" opt; do
    case $opt in
        e) e="$OPTARG";;
        b) b="$OPTARG";;
        n) n="$OPTARG";;
        s) s="$OPTARG";;
        i) i="$OPTARG";;
        \?) echo "Invalid option -$OPTARG" >&2; show_help; exit 1;;
    esac
done

# Remove the parsed options from the positional parameters
shift $((OPTIND-1))

# Check if the required taxid argument is provided
if [ $# -ne 1 ]; then
    echo "Invalid option" >&2; show_help; exit 1
fi

# Assign positional arguments to variables
taxid=$1

echo "Setting required variables:"
echo ""
echo "Selected taxid is: $taxid"
echo "Selected evalue threshold is: $e"
echo "Selected bitscore threshold is: $b"
echo "Selected max neighbours is: $n"
echo "Selected step is: $s"
echo "Selected number of iterations is: $i"

####---- DOWNLOADING FILES AND DATABASES ----####

echo ""
echo "Downloading required AlphaFold models:"
echo ""


# Define directory names
db_dir="foldseek_dbs"
prot_dir="proteomes"
results_dir="results"
taxid_dir="$results_dir/$taxid"

# Check if directories exist, create if necessary
mkdir -p "$db_dir" "$prot_dir" "$results_dir" "$taxid_dir" "$taxid_dir/SAFE"

# Check if files exist, download if necessary
if [ ! -f "$db_dir/afdb50sp" ]; then
    foldseek databases Alphafold/UniProt50-minimal "$db_dir/afdb50" tmp --remove-tmp-files 1
    foldseek databases Alphafold/Swiss-Prot "$db_dir/swissprot" tmp --remove-tmp-files 1

    # Merge the databases
    foldseek concatdbs "$db_dir/afdb50" "$db_dir/swissprot" "$db_dir/afdb50sp"
    foldseek concatdbs "$db_dir/afdb50_h" "$db_dir/swissprot_h" "$db_dir/afdb50sp_h"
    foldseek concatdbs "$db_dir/afdb50_ss" "$db_dir/swissprot_ss" "$db_dir/afdb50sp_ss"
    foldseek concatdbs "$db_dir/afdb50_ca" "$db_dir/swissprot_ca" "$db_dir/afdb50sp_ca"
else
    echo "Foldseek databases already downloaded"
fi

if [ -f "$taxid_dir/$taxid.m8" ] && [ -f "$taxid_dir/${taxid}_bh.m8" ] && [ -f "$taxid_dir/${taxid}_normalisation.m8" ] && [ -f "$taxid_dir/${taxid}_normalisation_bh.m8" ]; then
    echo "AlphaFold models of selected organism already downloaded and Foldseek results already generated"
elif [ ! -f "$prot_dir/$taxid.tar" ] && [ ! -f "$db_dir/afdb50sp$taxid" ]; then
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

    foldseek concatdbs "$db_dir/afdb50sp" "$db_dir/$taxid" "$db_dir/afdb50sp$taxid"
    foldseek concatdbs "$db_dir/afdb50sp_h" "$db_dir/$taxid_h" "$db_dir/afdb50sp${taxid}_h"
    foldseek concatdbs "$db_dir/afdb50sp_ss" "$db_dir/${taxid}_ss" "$db_dir/afdb50sp${taxid}_ss"
    foldseek concatdbs "$db_dir/afdb50sp_ca" "$db_dir/${taxid}_ca" "$db_dir/afdb50sp${taxid}_ca"
elif [ -f "$prot_dir/$taxid.tar" ] && [ ! -f "$db_dir/afdb50sp$taxid" ]; then
    foldseek createdb "$prot_dir/$taxid.tar" "$db_dir/$taxid" 

    foldseek concatdbs "$db_dir/afdb50sp" "$db_dir/$taxid" "$db_dir/afdb50sp$taxid"
    foldseek concatdbs "$db_dir/afdb50sp_h" "$db_dir/${taxid}_h" "$db_dir/afdb50sp${taxid}_h"
    foldseek concatdbs "$db_dir/afdb50sp_ss" "$db_dir/${taxid}_ss" "$db_dir/afdb50sp${taxid}_ss"
    foldseek concatdbs "$db_dir/afdb50sp_ca" "$db_dir/${taxid}_ca" "$db_dir/afdb50sp${taxid}_ca"
else
    echo "AlphaFold models of selected organism already downloaded"
fi

####---- RECIPROCAL BEST STRUCTURE HITS SEARCH ----####

echo ""
echo "Performing Reciprocal Best Structural Hits search:"
echo ""

if [ ! -f "$taxid_dir/$taxid.m8" ] || [ ! -f "$taxid_dir/${taxid}_bh.m8" ]; then
    foldseek easy-search --format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits" \
                         "$prot_dir/$taxid.tar" "$db_dir/afdb50sp$taxid" "$taxid_dir/$taxid.m8" tmp --threads 64

    foldseek easy-search --format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits" \
                         "$prot_dir/$taxid.tar" "$db_dir/$taxid" "$taxid_dir/${taxid}_normalisation.m8" tmp --threads 64 \
                         --exhaustive-search 1 --min-seq-id 0.9

    python3 get_besthits.py --input "$taxid_dir/$taxid.m8" --output "$taxid_dir/${taxid}_bh.txt" --eval_thr "$e" --bits_thr "$b"

    foldseek prefixid "$db_dir/afdb50sp${taxid}_h" "$db_dir/afdb50sp${taxid}.lookup" --tsv --threads 1
    awk 'NR == FNR {f[$1] = $1; next} $2 in f {print $1}' "$taxid_dir/${taxid}_bh.txt" "$db_dir/afdb50sp${taxid}.lookup" > "$db_dir/subset${taxid}.tsv"
    foldseek createsubdb "$db_dir/subset${taxid}.tsv" "$db_dir/afdb50sp$taxid" "$db_dir/subdb$taxid"
    foldseek createsubdb "$db_dir/subset${taxid}.tsv" "$db_dir/afdb50sp${taxid}_ss" "$db_dir/subdb${taxid}_ss"
    foldseek createsubdb "$db_dir/subset${taxid}.tsv" "$db_dir/afdb50sp${taxid}_ca" "$db_dir/subdb${taxid}_ca"
    rm "$db_dir/subset${taxid}.tsv"

    foldseek search "$db_dir/subdb$taxid" "$db_dir/afdb50sp$taxid" "$taxid_dir/${taxid}_bh" tmp -a 1 --threads 64
    foldseek convertalis --format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits" \
                         "$db_dir/subdb$taxid" "$db_dir/afdb50sp$taxid" "$taxid_dir/${taxid}_bh" "$taxid_dir/${taxid}_bh.m8"

    foldseek search "$db_dir/subdb$taxid" "$db_dir/subdb$taxid" "$taxid_dir/${taxid}_normalisation_bh" tmp -a 1 --threads 64 --exhaustive-search 1 --min-seq-id 0.9
    foldseek convertalis --format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits" \
                         "$db_dir/subdb$taxid" "$db_dir/subdb$taxid" "$taxid_dir/${taxid}_normalisation_bh" "$taxid_dir/${taxid}_normalisation_bh.m8"

    rm "$db_dir/subdb$taxid"*
    rm "$db_dir/afdb50sp$taxid"*
    rm "$db_dir/$taxid"*
    find "$taxid_dir" -type f ! -name "*.txt" ! -name "*.m8" -delete
else
    echo "Foldseek results already generated"
fi

psize=$(tar -tvf "$prot_dir/$taxid.tar" | grep ".gz" | wc -l)
python3 parse_m8.py --input "$taxid_dir/$taxid.m8" --input_normalisation "$taxid_dir/${taxid}_normalisation.m8" \
                    --input_bh "$taxid_dir/${taxid}_bh.m8" --input_normalisation_bh "$taxid_dir/${taxid}_normalisation_bh.m8" \
                    --proteome_size "$psize" --eval_thr "$e" --bits_thr "$b"

touch "$taxid_dir/${taxid}_nan.txt"

for j in $(seq 1 $i); do
    echo ""
    echo "Performing iteration $j"

    ####---- NETWORK GENERATION ----####
    echo ""
    echo "Creating RBSH network and identifying clusters of homologs:"
    echo ""

    neighbours=$((n + (s * (j - 1))))
    echo "Selected max number of neighbours for iteration $j is: $neighbours"
    echo ""
    python3 generate_network.py --input "$taxid_dir/$taxid.json" --input_bh "$taxid_dir/${taxid}_bh.json" --nan "$taxid_dir/${taxid}_nan.txt" --neighbours "$neighbours" \
                                --output "$taxid_dir/${taxid}_clusters_iter$j.txt" --edgelist "$taxid_dir/${taxid}_edgelist_iter$j.txt"

    ####---- ANNOTATION ----####
    echo ""
    echo "Annotating network's nodes"
    echo ""

    python3 retrieve_annotation.py --input "$taxid_dir/${taxid}_clusters_iter$j.txt" --output "$taxid_dir/${taxid}_annotation_iter$j.txt"

    ####---- SAFE ENRICHMENT AND STATISTICS COMPUTATION ----####
    echo "Performing SAFE analysis and computing statistics on new annotation"
    echo ""

    python3 SAFE_enrichment.py --input "$taxid_dir/${taxid}_annotation_iter$j.txt" --wd "$taxid_dir/" \
                               --edgelist "${taxid}_edgelist_iter$j.txt" --radius 0.5 --iteration $j --identifiers "$taxid_dir/${taxid}_normalisation.m8" \
                               --output "$taxid_dir/${taxid}_NEWannotation" --outfig "$taxid_dir/${taxid}_barcharts" --nan "$taxid_dir/${taxid}_nan.txt"

    if [ ! -s "$taxid_dir/${taxid}_nan.txt" ]; then
        echo "All nan2nan proteins have been annotated. Stopping at iteration $j."
        break
    else
        line_count=$(wc -l < "$taxid_dir/${taxid}_nan.txt")
        if [ "$j" -lt "$i" ]; then
            echo "Found $line_count nan2nan IDs in total in the target proteome... proceeding to next iteration."
        else
            echo "Found $line_count nan2nan IDs in total in the target proteome."
        fi
    fi
done

rm -r "$taxid_dir/SAFE"

echo ""
echo "WASP run completed!"
