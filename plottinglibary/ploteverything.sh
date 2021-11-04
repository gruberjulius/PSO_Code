#!/bin/bash

# to be called from automation folder
# if moved to different location path has to be changed
for filename in helperscripts/automation/*.txt; do
    python3 plottingwithclasses.py < ${filename}
done