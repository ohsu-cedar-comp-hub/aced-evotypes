#!/usr/bin/env python

import json
import argparse
import re
import random
import os
from minio import Minio

def get_args():
    parser = argparse.ArgumentParser(description="Script to generate json files for each sample")
    parser.add_argument("-m", "--manifest", help="absolute path to manifest file with matched tumor/normal samples", type=str, required=True)
    parser.add_argument("-g", "--global_params", help="absolute path to global_params.json", type=str, required=True)
    parser.add_argument("-o", "--output", help="absolute path to output directory where json files will be stored", type=str, required=True)
    return parser.parse_args()

args = get_args()
manifest = args.manifest
params = args.global_params
output_dir = args.output

# Initialize MinIO client
client = Minio(
    "rgw.ohsu.edu",
    access_key="<public key>",
    secret_key="<private key>",
    secure=True
)

bucket="gdc-esca"
# Function to find file location in the bucket
def find_file_location(file_name):
    objects = client.list_objects(bucket, recursive=True)
    for obj in objects:
        if file_name in obj.object_name:
            return obj.object_name  # Return the first match found
    return None  # Return None if not found

# Load the existing params in global_params.json
with open (params, 'r') as file:
    data = json.load(file)

# Random number parameters
lower_bound=1
upper_bound=1000
unique_numbers = set() # keep track of numbers used

# Getting sample information from the manifest file
with open(manifest, 'r') as file: 
    for line in file:
        line = line.strip().split()

        if line == "":
            break
        else:
            # col1 file
            file_1 = line[0]
            file_1_basename = os.path.splitext(file_1)[0]
            file_1_path = find_file_location(file_1)
            
            # col2 file
            file_2 = line[1]
            file_2_basename = os.path.splitext(file_2)[0]
            file_2_path = find_file_location(file_2)

            # getting Case ID
            x = re.findall(f'DO[0-9]+', line[0])
            case_id = x[0]

            # generate random number
            while True:
                number_var = random.randint(lower_bound, upper_bound)
                if number_var not in unique_numbers:
                    unique_numbers.add(number_var)
                    break

            # generate random number
            while True:
                number_cgpwgs = random.randint(lower_bound, upper_bound)
                if number_cgpwgs not in unique_numbers:
                    unique_numbers.add(number_cgpwgs)
                    break

            # adding local params to global params
            data["number_var"]=str(number_var)
            data["number_cgpwgs"]=str(number_cgpwgs)
            data["case_id"]=case_id
            data["file_1"]=file_1
            data["file_1_basename"]=file_1_basename
            data["file_1_path"]=f's3://gdc-esca/{file_1_path}'
            data["file_2"]=file_2
            data["file_2_basename"]=file_2_basename
            data["file_2_path"]=f's3://gdc-esca/{file_2_path}'

            writing JSON to a file
            with open(f'{output_dir}/{case_id}.json', 'w') as out:
                json.dump(data, out, indent=4)
            out.close()
            print(f'JSON data has been stored in "{case_id}.json"')
file.close()
