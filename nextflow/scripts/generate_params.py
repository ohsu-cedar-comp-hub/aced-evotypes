#!/usr/bin/env python

import json
import argparse
import re
import random
from minio import Minio

def get_args():
    parser = argparse.ArgumentParser(description="Script to generate a json file per sample")
    parser.add_argument("-m", "--manifest", help="Absolute path to the file with matched tumor/normal samples", type=str, required=True)
    parser.add_argument("-g", "--global_params", help="Absolute path to global_params.json", type=str, required=True)
    parser.add_argument("-o", "--output", help="Absolute path to a directory where output json files will go", type=str, required=True)
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
            return obj.object_name
    return None # if file not found

# Load the existing params in global_params.json
with open (params, 'r') as file:
    data = json.load(file)

# Random number parameters (random number will be used when creating scratch dir before each run)
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
            # normal file
            normal_file = line[0]
            normal_basename = line[0].rsplit('.', 1)[0]
            normal_path = find_file_location(normal_file)
            
            # tumor file
            tumor_file = line[1]
            tumor_basename = line[1].rsplit('.',1)[0]
            tumor_path = find_file_location(tumor_file)

            # getting Case ID
            x = re.findall(f'DO[0-9]+', line[0])
            case_id = x[0]

            # generate random number
            while True:
                random_number = random.randint(lower_bound, upper_bound)
                if random_number not in unique_numbers:
                    unique_numbers.add(random_number)
                    break

            # adding local params to global params
            data["number"]=str(random_number)
            data["case_id"]=case_id
            data["normal_file"]=normal_file
            data["normal_basename"]=normal_basename
            data["normal_path"]=f's3://gdc-esca/{normal_path}'
            data["tumor_file"]=tumor_file
            data["tumor_basename"]=tumor_basename
            data["tumor_path"]=f's3://gdc-esca/{tumor_path}'

            # writing JSON to a file
            with open(f'{output_dir}/{case_id}.json', 'w') as out:
                json.dump(data, out, indent=4)
            out.close()
            print(f'JSON data has been stored in "{case_id}.json"')
            samples = {}
file.close()