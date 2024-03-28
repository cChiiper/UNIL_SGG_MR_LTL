#!/bin/bash

# Define inputs chr6:25'000'000–37'000'000
chromosome=6 
start_pos=25000000 # Start position
end_pos=37000000 # End position
input_file="/uk10k.autosomal.bim"
output_file="/HLA_uk10k_rsids.txt"

# Use awk to process the file
awk -v chr="$chromosome" -v start="$start_pos" -v end="$end_pos" '{
    split($2, rsids, ";"); # Split the rsid field on semicolon
    if ($1 == chr && $4 >= start && $4 <= end) {
        for (id in rsids) {
            print rsids[id]
        }
    }
}' "$input_file" > "$output_file"