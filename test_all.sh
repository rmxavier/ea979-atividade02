#!/bin/bash

for file in test_all_commands test_defaults test_geometry test_transforms; do
    python3 draw_3d_model.py "${file}.dat" "${file}.ppm"

    if [ -s "${file}.ppm" ]; then
        diff -q "${file}.ppm" "${file}-ref.ppm"
        if [ "$?" == "0" ]; then
            rm -f "${file}-diff.ppm"
        else
            # In systems that have no pamarith, change for:
            # pnmarith -difference "${file}.ppm" "${file}-ref.ppm" > "${file}-diff.ppm"
            pamarith -xor "${file}.ppm" "${file}-ref.ppm" > "${file}-diff.ppm"
        fi
    else
        echo "error creating ${file}.ppm"
    fi
done
