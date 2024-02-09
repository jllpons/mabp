#!/usr/bin/env bash

# Function to get PDB code either from command-line argument or stdin
get_pdb_code() {
    if [ "$#" -eq 1 ]; then
        echo "$1"
    else
        read pdb_code
        echo "$pdb_code"
    fi
}

# Attempt to get PDB code
pdb_code=$(get_pdb_code "$@")

# Check if PDB code was provided
if [ -z "$pdb_code" ]; then
    echo "Usage: $0 [PDB_CODE] or echo PDB_CODE | $0" >&2
    exit 1
fi

# Start the mapping job
jobId=$(curl -s --form 'from="PDB"' \
                  --form 'to="UniProtKB"' \
                  --form 'ids="'$pdb_code'"' \
                  https://rest.uniprot.org/idmapping/run | jq -r '.jobId') >&2

if [ -z "$jobId" ]; then
    echo "Failed to submit job." >&2
    exit 1
fi

echo "Job ID: $jobId" >&2

# Initialize job status
status="RUNNING"
max_attempts=30
attempts=0

# Wait for the job to complete
while [ "$status" == "RUNNING" ] && [ $attempts -lt $max_attempts ]; do
    echo "Checking job status..." >&2
    status_response=$(curl -s https://rest.uniprot.org/idmapping/status/$jobId)
    echo "Raw status response: $status_response" >&2  # Debug: Print raw JSON response
    status=$(echo "$status_response" | jq -r '.jobStatus')
    echo "Parsed status: $status" >&2  # Debug: Print parsed status
    if [ "$status" == "FINISHED" ]; then
        break
    elif [ "$status" == "FAILURE" ]; then
        echo "Job failed." >&2
        exit 2
    fi
    sleep 10 # Wait for 10 seconds before checking again
    ((attempts++))
done

if [ $attempts -eq $max_attempts ]; then
    echo "Max attempts reached without completing." >&2
    exit 3
fi

# Check final status and retrieve the results if successful
if [ "$status" == "FINISHED" ]; then
    echo "Retrieving results..." >&2
    results=$(curl -s https://rest.uniprot.org/idmapping/results/$jobId)
    echo "$results" | jq -r '.results[].to'  # Adjust based on the actual JSON structure
else
    echo "Job did not finish successfully." >&2
    exit 4
fi

