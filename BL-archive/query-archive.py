#%%
"""
Code Purpose: Cross reference Breakthrough Listen Observations with Known Radio Flare Stars 
Author: Owen A. Johnson 
Last Updated: 2024-04-24 
"""
import requests
from tqdm import tqdm
import pandas as pd
import time

TARGETS_URL = "http://seti.berkeley.edu/opendata/api/list-targets"

response = requests.get(TARGETS_URL)

if response.status_code == 200:
    trgt_list = response.json()
    print(f"Number of targets on BL API: {len(trgt_list)}")
else:
    print(f"Error: {response.status_code}")
    trgt_list = []

API_URL = "http://seti.berkeley.edu/opendata/api/query-files"
all_entries = []
counter = 0

for target in tqdm(trgt_list[1:]):
    payload = {
        "target": target,
        "limit": 20000
    }
    
    response = requests.get(API_URL, params=payload)
    if response.status_code == 200:
        query_data = response.json()
        metadata = query_data['data']
        for entry in metadata:
            extracted = {
                "telescope": entry.get("telescope", "unknown"),
                "center_freq": entry.get("center_freq", "unknown"),
                "ra_deg": entry.get("ra", "unknown"),
                "dec_deg": entry.get("decl", "unknown"),
                "data_url": entry.get("url", "unknown"),
                "utc": entry.get("utc", "unknown"),
                "mjd": entry.get("mjd", "unknown"),
                "source_name": entry.get("target", "unknown")
            }
            all_entries.append(extracted)
    else:
        print(f"Error: {response.status_code} for target {target}")
        
    counter += 1
    if counter % 100 == 0:
        df = pd.DataFrame(all_entries)
        df.to_csv("BL_opendata_metadata.csv", index=False)
    
    time.sleep(0.5)

df = pd.DataFrame(all_entries)
df.to_csv("BL_opendata_metadata.csv", index=False)
