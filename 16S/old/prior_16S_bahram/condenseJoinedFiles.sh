#!/bin/bash

for dir in *; do
    if test -d "$dir"; then
        (
            cd $dir
            for file in *; do
                echo $file
                if [[ "$file" == *.join.fastq ]]; then
                    newfile=$dir.${file##*.}
                    cp "$file" "../$newfile"
                fi
            done
        )
    fi
done
