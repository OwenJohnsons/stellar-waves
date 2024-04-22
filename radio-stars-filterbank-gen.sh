source="$1" # e.g. "J0835-4510"
path="${2:-../}" # e.g. "/mnt/data/lofar/udp/2021-08-01"

zst_files=(${path}/*.zst)
start_time=$(echo "${zst_files[0]}" | awk -F '[/_T]' '{print $9}' | sed 's/..$//')
date=$(echo "${zst_files[0]}" | cut -d'/' -f8 | cut -d'.' -f3 | cut -dT -f1)
time=$(echo "${zst_files[0]}" | cut -d'/' -f8 | cut -dT -f2 | cut -d'.' -f1)

echo " --- Processing filterbanks for $source ---"
echo " Path: $path"
echo " Number of .zst files: $(ls -1a "${path}"/udp*.zst | wc -l)"
echo " Obs Date: $date"
echo " Obs Time: $time"

pattern="udp_[0-9]*\.ucc1"
matched_part=$(echo "$zst_files[0]" | grep -oE "$pattern")

if [[ -n "$matched_part" ]]; then
    # Replace the matched part with 'udp_%d.ucc1'
    modified_string="${zst_files[0]/$matched_part/udp_1613%d.ucc1}"
else
    echo "No matching pattern found in the string."
fi

# lofar_udp_extractor -i ${modified_string} -u 4 -p 10 -o "${path}/$source"_"$start_time".dada | tee -a filgen_output_"$source"_"$start_time".log
lofar_udp_extractor -i ${modified_string} -u 4 -p 10 -o "/mnt/ucc4_data1/data/Owen/dada_dump/$source"_"$start_time".dada | tee -a filgen_output_"$source"_"$start_time".log
