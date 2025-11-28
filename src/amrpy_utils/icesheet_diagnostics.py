####################################################################################

import os
import argparse
import subprocess
import pandas as pd

####################################################################################

# Set constants
ICE_DENSITY = 918.0
WATER_DENSITY = 1028.0
H_MIN = 1.0e-3

####################################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run BISICLES diagnostics on plot.*.hdf5 files.")
    parser.add_argument("--filetools_path", default= "./", help="Path to the filetools directory")
    parser.add_argument("--diagnostics_exec", default="diagnostics2d.Linux.64.CC.ftn.DEBUG.OPT.MPI.GNU.ex", help="Name of the diagnostics executable")
    parser.add_argument("--output_path", default = "./", help="Path to the output directory")
    parser.add_argument("--icesheet", choices=["AIS", "GrIS"], default="AIS", help='Select the icesheet: "AIS" for Antarctic Ice Sheet, "GrIS" for Greenland Ice Sheet.')
    parser.add_argument("--directory", help="Directory containing plot.*.hdf5 files")
    parser.add_argument("--csv_file", help="Name of the output CSV file")
    parser.add_argument("--mask", nargs="?", default=0, help="Optional mask file (default: 0 for no mask)")
    parser.add_argument("--mask_no_start", type=int, default=0, help="ID of the first subdomain to calculate - 0 is the entire domain (default: 0)")
    parser.add_argument("--mask_no_end", type=int, default=16, help="ID of the last subdomain to calculate (default: 16)")
    return parser.parse_args()

def validate_directory(directory):
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory does not exist: {directory}")
    return os.path.abspath(directory)

def construct_command(plot_file, csv_file, append=False, mask_file=None, mask_no_start=0, mask_no_end=16):
    cmd = f"{os.path.join(FILETOOLS_PATH, DIAGNOSTICS_EXEC)} plot_file={plot_file} out_file={os.path.join(OUTPUT_PATH, csv_file)}"
    if append:
        cmd += " -append"
    cmd += f" ice_density={ICE_DENSITY} water_density={WATER_DENSITY} h_min={H_MIN}"
    if mask_file:
        cmd += f" mask_file={mask_file} mask_no_start={mask_no_start} mask_no_end={mask_no_end}"
    return cmd

####################################################################################

args = parse_arguments()

# Validate and process directory
icesheet = args.icesheet
directory = validate_directory(args.directory)
csv_file = args.csv_file
mask_file = None if args.mask == "0" else args.mask
masked = mask_file is not None
mask_no_start = args.mask_no_start
mask_no_end = args.mask_no_end

print(f"Running icesheet diagnostics in directory: {directory}")
if masked:
    print(f"Using mask file: {mask_file}")

# Check to see if there is already a csv file for this ice sheet, suite, mask combnination
csv_exists = os.path.exists(os.path.join(OUTPUT_PATH, csv_file))

# If there is, then only process the files that have not yet been processed, and append to the csv file
if csv_exists:
    csv_data = pd.read_csv(os.path.join(OUTPUT_PATH, csv_file))
    processed_files = csv_data["filename"].values.tolist()
    processed_files = list(dict.fromkeys(processed_files)) # Remove duplicates
else:
    processed_files = []

all_files = sorted(os.listdir(directory))
all_files = [file for file in all_files if file.startswith("plot.")]
all_files_abs = [os.path.join(directory, file) for file in all_files]

# Remove files that have already been processed
if len(processed_files) > 0:
    print(f"{len(processed_files)} files have already been processed. Removing from list of files to process.")
    for file in processed_files:
        all_files_abs.remove(file)

# If all files have been processed then exit
if len(all_files_abs) == 0:
    print("All relevant hdf5 files have already been processed. Exiting now.")

else:
    print(f"No. of files left to process is: {len(all_files_abs) - len(processed_files)}")

    # Process files
    count = 0

    for infile in all_files_abs:

        plot_file = infile
            
        append = csv_exists or count > 0
        
        diags_command = construct_command(plot_file, csv_file, append=append, mask_file=mask_file, mask_no_start=mask_no_start, mask_no_end=mask_no_end)

        print(f"Running diagnostics for file: {infile}")

        subprocess.run(diags_command, shell=True)

        count += 1

    # Summary
    if count == 0:
        print("No BISICLES *.hdf5 files found needing processing in the directory.")
    else:
        print(f"{count} *.hdf5 files processed. Diagnostics written to {csv_file}.")