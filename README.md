# HapSim-Lite

Proof-of-concept lightweight Markov-Chain N-ploidy variant simulation with preserved MAF and LD. 

## Overview

Simulation of variant data is a key component of method development and validation. 
Many methods have been developed and released in the past 15 years, and mature tools like HapGen have provided robust methods for simulating small to biobank scale variation data with realistic linkage disequilibrium and minor allele frequencies. 

HapSim-Lite is intended to serve as a lightweight option for prototyping and production use with a reasonable amount of tuning. 
Its CLI outputs to VCF natively, removing complex dependencies and provides straightforward simulation of arbitrary numbers of samples and ploidy with phased or unphased output, even if the input LD and MAF is calculated from diploid or haploid data. 

## Installation

Install from source: 

```bash
pip install https://github.com/sadams2013/hapsim-lite#egg=hapsim_lite
```

# Usage

## CLI

### Generating input frequencies and LD data from a reference panel (e.g. 1000 genomes)

HapSim-Lite has a CLI (`hapsim-lite`) that is intended to work with frequency and LD outputs from [Plink2](https://www.cog-genomics.org/plink/2.0/) 

Chromosome and whole-genome pvar/psam/pgen files can be obtained here: https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg or generated from your own VCF files. Note that you should split multiallelic markers before generating. 

**Plink input files MUST use variant IDs formatted as chr_pos_ref_alt.** This can be set with (or adapted for whatever the input is (e.g. VCF))

```bash

plink/plink2 \
    --pfile base/plink/prefix \
    --set-all-var-ids @_#_\$r_\$a \
    --new-id-max-allele-len 1000 \
    --out preprocessed \
    --make-pgen
```

The following plink commands will generate required input files for hapsim-lite. The maf and hwe filters are optional, but will help with runtime by removing rare and/or low-quality variants that are less likely to perform well in the simulated output. They may be adjusted or omitted depending on your objective.

```bash

# minor allele frequencies

plink2 \
    --pfile preprocessed \
    --freq \
    --maf 0.05 \
    --hwe 1e-6 0.001 \
    --out FREQ

## Output file: FREQ.afreq

# linkage disequilibrium

plink2 \
    --pfile preprocessed \
    --maf 0.05 \
    --hwe 1e-6 0.001 \
    --r-unphased \
    --ld-window-kb 20 \
    --ld-window-r2 0.1 \
    --out LD

## Output file: LD.vcor

```

### Running the CLI

Full usage available with `hapsim-lite -h`

A suitable starting run to simulate 100 phased diploid samples might look like: 

```bash
hapsim-lite \
    -f FREQ.afreq \
    -l LD.vcor \
    -n 100 > simulated.vcf
```

hapsim-lite CLI outputs directly to stdout, which should be piped through `bcftools view` or `bgzip`. 

# TODO

- Examples and documentation for interaction with the population and simulation API
- Allow user-specified seed for deterministic outputs
- Parse variant info from pvar to relax the variant ID format requirement
- Further hyperparameter tuning to improve recovery of LD signal

