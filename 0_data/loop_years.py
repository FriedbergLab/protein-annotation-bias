import os
import re
import subprocess
import sys

def get_filtered_files(directory, pattern):
    files = os.listdir(directory)
    filtered_files = sorted([f for f in files if pattern.match(f)])
    return sorted(filtered_files, key=lambda x: int(pattern.match(x).group(1)))


def create_file_paths(directory, files):
    return [os.path.join(directory, f) for f in files[:10]]


def run_propagation(gaf_file, year, pivot_year, pivot_obo_file):
    args = [
        'python3', '0_data/extract_GOA.py',
        str(year), gaf_file, str(pivot_year), pivot_obo_file, 'EXPEC', '1'
    ]
    print(f'Running command: {" ".join(args)}')
    subprocess.run(args)


def main():
    if len(sys.argv) != 3:
        print('Usage: python3 loop_years.py <gaf_dir> <obo_dir>')
        sys.exit(1)
        
    gaf_dir, obo_dir = sys.argv[1], sys.argv[2]
    gaf_pattern = re.compile(r'.*([0-9]{2})goa_human\.gaf$')
    obo_pattern = re.compile(r'([0-9]{2})go\.obo$')

    filtered_gaf_files = get_filtered_files(gaf_dir, gaf_pattern)
    filtered_obo_files = get_filtered_files(obo_dir, obo_pattern)

    gaf_file_paths = create_file_paths(gaf_dir, filtered_gaf_files)
    obo_file_paths = create_file_paths(obo_dir, filtered_obo_files) # sorted based on increasing year

    years = list(range(2022, 2012, -1))
    pivot_year = years[0]
    pivot_obo_file = obo_file_paths[-1]

    # Run 2022 first
    run_propagation(gaf_file_paths[-1], years[0], pivot_year, pivot_obo_file)

    # Then run the rest of the years in descending order
    for gaf_file, year in zip(gaf_file_paths[-2::-1], years[1:]):
        run_propagation(gaf_file, year, pivot_year, pivot_obo_file)


if __name__ == "__main__":
    main()
