source_="$1" # e.g. "CR_Draconis"
source=$(echo $source_ | sed 's/_/ /g')
path="${2:-../}" # e.g. "/mnt/data/lofar/udp/2021-08-01"

zst_files=(${path}/*.zst)
start_time=$(echo "${zst_files[0]}" | awk -F '[/_T]' '{print $9}' | sed 's/..$//')
date=$(echo "${zst_files[0]}" | cut -d'/' -f8 | cut -d'.' -f3 | cut -dT -f1)
time=$(echo "${zst_files[0]}" | cut -d'/' -f8 | cut -dT -f2 | cut -d'.' -f1)
formatted_date=$(echo $start_time | sed 's/\(....\)\(..\)\(..\)\(..\)\(..\)/\1-\2-\3T\4:\5:/; s/\(.*\):/\1\:/')
RA=$(python -c "from astropy.coordinates import SkyCoord; c = SkyCoord.from_name('$source'); print(c.ra.to_string(unit='hourangle', sep=':', precision=2))")
DEC=$(python -c "from astropy.coordinates import SkyCoord; c = SkyCoord.from_name('$source'); print(c.dec.to_string(unit='degree', sep=':', precision=2))")
tstartmjd=$(python -c "from astropy.time import Time; print(Time('$formatted_date').mjd)" )


echo " =========== Processing filterbanks for $source =========== "
echo " Path: $path"
echo " Number of .zst files: $(ls -1a "${path}"/udp*.zst | wc -l)"
echo " Obs Date: $date"
echo " Obs Time: $time"
echo " RA: $RA"
echo " DEC: $DEC"
echo " MJD Start: $tstartmjd"
echo " =========================================================== "

datetime="${obs_date} ${obs_time}"

pattern="udp_[0-9]*\.ucc1"
matched_part=$(echo "$zst_files[0]" | grep -oE "$pattern")

if [[ -n "$matched_part" ]]; then
    # Replace the matched part with 'udp_%d.ucc1'
    modified_string="${zst_files[0]/$matched_part/udp_1613%d.ucc1}"
else
    echo "No matching pattern found in the string."
fi

lofar_udp_extractor \
    -i "${modified_string}" \
    -u 4 \
    -p 154 \
    -a "-fch1 200 -fo -0.1953125 -source ${source_} -ra ${RA} -dec ${DEC}" \
    -o "/mnt/ucc4_data2/data/Owen/raw-radio-stars/${source_}_${formatted_date}_raw_S%d.fil" | tee -a "./logs/filgen_output_${source_}_${formatted_date}.log"

# # --- Fil Generation ---
for i in {0..3}; do
    if [ ! -f "/mnt/ucc4_data2/data/Owen/raw-radio-stars/${source_}_${formatted_date}_S${i}.fil" ]; then
        digifil -b-32 -t 128 -c -I 0 "/mnt/ucc4_data2/data/Owen/raw-radio-stars/${source_}_${formatted_date}_raw_S${i}.fil" \
        -o "/mnt/ucc4_data2/data/Owen/raw-radio-stars/${source_}_${formatted_date}_S${i}.fil"
        chmod 777 "/mnt/ucc4_data2/data/Owen/raw-radio-stars/${source_}_${formatted_date}_S${i}.fil"
    else
        echo "File already exists: /mnt/ucc4_data2/data/Owen/raw-radio-stars/${source_}_${formatted_date}_S${i}.fil skipping digifil"
    fi
done

# --- Clean Up ---
rm -f /mnt/ucc4_data2/data/Owen/raw-radio-stars/${source}_${formatted_date}_raw_S*.fil
