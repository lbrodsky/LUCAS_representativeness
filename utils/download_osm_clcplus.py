import os
import sys
import requests
import argparse
import py7zr
from pathlib import Path

from utils import country_codes

def download_file(url, dst_dir):
    output_file = Path(dst_dir, url.split('/')[-1])

    response = requests.get(url)
    with open(output_file, "wb") as f:
        f.write(response.content)

    return output_file

def extract_data(archive_path, dst_dir):
    with py7zr.SevenZipFile(archive_path, mode='r') as archiv:
        archiv.extractall(path=dst_dir)

def main(dst_dir, year):
    if year == 2018:
        zenodo_id = 15648881
    elif year == 2022:
        zenodo_id = 15639393
    else:
        sys.exit(f"ERROR: Year {year} not supported")

    if not dst_dir.exists():
        dst_dir.mkdir(parents=True, exist_ok=True)
        
    for cntr in sorted(country_codes(year).keys()):
        print(f"Downloading {cntr}...", file=sys.stderr)
        data_file = download_file(f"https://zenodo.org/records/{zenodo_id}/files/{cntr}_osm_clcplus_{year}.7z",
                                  dst_dir)
        extract_data(data_file, Path(dst_dir) / str(year))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dst_dir', help='Destination directory.')
    parser.add_argument('--year', help='Year to be downloaded.', type=int, default=2018)

    args = parser.parse_args()

    main(Path(args.dst_dir), args.year)
