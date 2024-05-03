source_="$1" # e.g. "CR_Draconis"
source=$(echo $source_ | sed 's/_/ /g')
path="${2:-../}" # e.g. "/mnt/data/lofar/udp/2021-08-01"

zst_files=(${path}/*.zst)
start_time=$(echo "${zst_files[0]}" | awk -F '[/_T]' '{print $9}' | sed 's/..$//')
date=$(echo "${zst_files[0]}" | cut -d'/' -f8 | cut -d'.' -f3 | cut -dT -f1)
time=$(echo "${zst_files[0]}" | cut -d'/' -f8 | cut -dT -f2 | cut -d'.' -f1)
formatted_date=$(echo $start_time | sed 's/\(....\)\(..\)\(..\)\(..\)\(..\)/\1:\2:\3_\4:\5:/; s/\(.*\):/\1\./')

echo " =========== Processing filterbanks for $source =========== "
echo " Path: $path"
echo " Number of .zst files: $(ls -1a "${path}"/udp*.zst | wc -l)"
echo " Obs Date: $date"
echo " Obs Time: $time"

datetime="${obs_date} ${obs_time}"

pattern="udp_[0-9]*\.ucc1"
matched_part=$(echo "$zst_files[0]" | grep -oE "$pattern")

if [[ -n "$matched_part" ]]; then
    # Replace the matched part with 'udp_%d.ucc1'
    modified_string="${zst_files[0]/$matched_part/udp_1613%d.ucc1}"
else
    echo "No matching pattern found in the string."
fi

lofar_udp_extractor -i ${modified_string} -u 4 -p 154 -o "/mnt/ucc4_data2/data/Owen/raw-radio-stars/$source_"_"$formatted_date"_S%d.dada | tee -a ./logs/filgen_output_"$source_"_"$formatted_date".log
grep "Start time" ./logs/filgen_output_"$source_"_"$formatted_date".log | awk '{print $3}' > ./logs/start_time_"$source_"_"$formatted_date".log
tstart_true=$(cat ./logs/start_time_"$source_"_"$formatted_date".log | head -n 1)

# 4 pol, 488 beams, 195312.5 samples per second -> 23828125 bytes per second due to 16x time scrunch 
bytes=$(echo $tstart_true | awk -F. '{print "0."$2}' | awk '{printf("%d\n", $1 * 23828125) }')
echo " Bytes: $bytes"

# # --- Header Manipulation ---
for i in {0..3}; do
    hdr_loc="/mnt/ucc4_data2/data/Owen/raw-radio-stars/${source_}_${formatted_date}_S${i}.hdr"
    cp ./hdrtemplate.hdr "$hdr_loc"
    date=$(echo "$tstart_true" | awk -F'T' '{print $1}')
    time=$(echo "$tstart_true" | awk -F'T' '{print $2}' | awk -F. '{print $1}')
    tstartmjd=$(echo "$tstart_true" | python -c "from astropy.time import Time; import sys; time = sys.stdin.readline().strip(); print(Time(time).mjd)" )
    RA=$(python -c "from astropy.coordinates import SkyCoord; c = SkyCoord.from_name('$source'); print(c.ra.to_string(unit='hourangle', sep=':', precision=2))")
    DEC=$(python -c "from astropy.coordinates import SkyCoord; c = SkyCoord.from_name('$source'); print(c.dec.to_string(unit='degree', sep=':', precision=2))")
    sed -i "s|isotstart|${date}-${time}|g" "$hdr_loc"
    sed -i "s|isooffset|${bytes}|g" "$hdr_loc"
    sed -i "s|src_name|${source}|g" "$hdr_loc"
    sed -i "s|racoords|${RA}|g" "$hdr_loc"
    sed -i "s|decoords|${DEC}|g" "$hdr_loc"
done


# --- Fil Generation ---
for i in {0..3}; do
    digifil -b 8 -t 128 "/mnt/ucc4_data2/data/Owen/raw-radio-stars/${source}_${formatted_date}_S${i}.dada" -o "/mnt/ucc4_data2/data/Owen/raw-radio-stars/${source}_${formatted_date}_S${i}.fil"
done
