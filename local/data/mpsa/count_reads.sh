#!/bin/bash

# The pattern of the files to match
pattern='results_files/*.txt'

for file in *$pattern*
do
    if [ -f "$file" ]; then
        sum_inc=$(awk -F '\t' 'NR > 1 {sum += $6} END {print sum}' "$file")
        sum_tot=$(awk -F '\t' 'NR > 1 {sum += $7} END {print sum}' "$file")
        IFS='.' read -ra str_array <<< "$file"
        name="${str_array[1]}"
        echo -e "$name\tinc\t$sum_inc"
        echo -e "$name\ttot\t$sum_tot"
    fi
done

# The pattern of the files to match
pattern='cipher_files/*.txt'

for file in *$pattern*
do
    if [ -f "$file" ]; then
        sum=$(awk -F '\t' 'NR > 1 {sum += $3+$4} END {print sum}' "$file")
        IFS='.' read -ra str_array <<< "$file"
        name="${str_array[1]}"
        echo -e "$name\tssbc\t$sum"
    fi
done


