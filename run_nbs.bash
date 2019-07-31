#!/bin/bash

set -e

RESULTSDIR="results/notebooks"

mkdir -p $RESULTSDIR

declare -a nbs=("wt_neut_and_bind.ipynb"
                "analyze_map.ipynb")

for nb in "${nbs[@]}"
do
    echo "Running $nb"

    jupyter nbconvert \
        --to notebook \
        --execute \
        --inplace \
        --ExecutePreprocessor.timeout=-1 \
        $nb

    echo "Converting $nb to Markdown"
    jupyter nbconvert \
        --output-dir $RESULTSDIR \
        --to markdown \
        $nb
done
