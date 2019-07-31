# Mutational antigenic profiling of Perth/2009 H3 HA against Vietnam cohort sera

## Overview
This repository contains mutational antigenic profiling of serum samples collected in Vietnam as described in the following paper:
[Structure of general-population antibody titer distributions to influenza A virus](https://www.nature.com/articles/s41598-017-06177-0).
In that study, they measured binding titers to the HA1 domains of a 2009 pdmH1N1 virus and a 2009 and 2011 H3N2 virus.
Here we perform mutational antigenic profiling of some of these samples against the HA of the Perth/2009 H3N2 strain.

The sample collection and binding titer measurements were led by [Maciek Boni](https://bio.psu.edu/directory/mfb9) and numerous collaborators in Vietnam (see [here](https://www.nature.com/articles/s41598-017-06177-0)).
The mutational antigenic profiling and analysis was done by Rachel Eguia, Juhye Lee, and [Jesse Bloom](https://research.fhcrc.org/bloom/en.html).

## Quick overview of results
Look at the Markdown output of the Jupyter notebooks for a quick overview of the results:

  - [results/notebooks/wt_neut_and_bind.md](results/notebooks/wt_neut_and_bind.md) has neutralization assays of the sera versus the wildtype Perth/2009 HA, and comparison to the binding titers.
  
  - [results/notebooks/analyze_map.md](results/notebooks/analyze_map.md) has the results of the mutational antigenic profiling.

## Configuring and running the analysis
The configuration for the analysis is in [config.yaml](config.yaml), which in turns points to various other lists of samples and input data.
This configuration file is self-explanatory.

The analysis is performed by a series of Jupyter notebooks:

  - [wt_neut_and_bind.ipynb](wt_neut_and_bind.ipynb)

  - [analyze_map.ipynb](analyze_map.ipynb)

To execute these notebooks and generate the Markdown output referred to in the subsection above, runs the script [run_nb.bash](run_nbs.bash) with:

    ./run_nbs.bash

## Input data
All input data is in the [./data/](data) subdirectory.
Specifically:

  - [data/serum_info.csv](data/serum_info.csv): information on the serum samples; columns are:

    - *serum*: name of serum sample

    - *serum_group*: group of samples to which serum belongs

    - *age*: age of individual from whom serum was obtained

    - *H1_2009_binding*: binding titer to A/California/6/2009 (pdmH1N1) HA

    - *H3_2009_binding*: binding titer to A/Victoria/210/2009 (H3N2) HA

    - *H3_2011_binding*: binding titer to A/Victoria/361/2011 (H3N2) HA

    - *serum_sample_subset*: subset group to which serum belongs in initial partitioning used by Maciek to select samples for this study

  - [data/wt_neut_config](data/wt_neut_config.yaml): information on plate reader data for the neutralization assays of the wildtype Perth/2009 HA against the sera.

## Results
The results are placed in the [./results/](results) subdirectory.
Only some of these results are tracked in this GitHub repo since some of the files are quite large.

