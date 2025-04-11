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

def main(dst_dir):
    for cntr in sorted(country_codes.keys()):
        print(f"Downloading {cntr}...", file=sys.stderr)
        data_file = download_file(f"https://zenodo.org/records/15039461/files/{cntr}_osm_clcplus_2018.7z",
                                  dst_dir)
        extract_data(data_file, dst_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dst_dir', help='Destination directory.')

    args = parser.parse_args()

    if not Path(args.dst_dir).exists():
        os.makedirs(args.dst_dir)

    main(args.dst_dir)
