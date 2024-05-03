#!/bin/bash

mkdir ./output_fils/
chmod o+rw output_fils
chown -R 1000:1000 output_fils

lastpath="${1:-$(find /mnt/ucc1_recording2/data/observations/ -maxdepth 1 -type d -name "20*" | sort | tail -n 1)}"

mv ${lastpath}/*Sun357 ./

find ./20*Sun357/ -type f -name "udp_16130*zst" | while read -r portZero; do \
        formattedPortZero=${portZero/udp_16130/udp_1613\%d}; \
        date="$(basename "$portZero" | sed -re 's/.*([0-9-]{10}T[0-9:]{8}).*/\1/g')"; \
        outputFilterbank=./output_fils/Sun357_"$date"_StokesVector\%d_D000.fil; \
	if [ -f "${outputFilterbank/\%d/0}.zst" ]; then continue; fi; \
	if [ -f "${outputFilterbank/\%d/0}" ]; then continue; fi; \
        lofar_udp_extractor \
                -p 164 \
                -a "-fch1 200 -fo -0.1953125 -source Sun357" \
                -d "0,0,SUN" \
                -c "LBA,54:453:2 HBA,54:453:2 HBA,54:229:2" \
                -i "$formattedPortZero" \
                -o "${outputFilterbank/_D000.fil/.fil}" | tee -a ./output_fils/Sun_"$date"_processing.log;
	digifil -c -b-32 -I 0 -t 16 "${outputFilterbank/\%d_D000.fil/0.fil}" -o "${outputFilterbank/\%d/0}"
	digifil -c -b-32 -I 0 -t 16 "${outputFilterbank/\%d_D000.fil/1.fil}" -o "${outputFilterbank/\%d/1}"
	rm "${outputFilterbank/\%d_D000/0}" "${outputFilterbank/\%d_D000/1}"
	chmod o+r ./*/*
done

echo "Processing complete; compressing data."
for fil in ./output_fils/*"$date"*_D000.fil; do \
        zstd -3 -T8 --rm "$fil" -o "$fil".zst; \
done

chown -R 1000:1000 .
