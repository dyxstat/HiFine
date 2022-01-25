# HiFine: integrating Hi-c-based and shotgun-based methods to reFine binning of metagenomic contigs

## Introduction
HiCBin is a novel binning pipeline to refine the binning results of metagenomic contigs by integrating both Hi-C-based and shotgun-based binning tools. HiFine designs a strategy of fragmentation for the original bin sets derived from the Hi-C-based and shotgun-based binning methods, which considerably increases the purity of initial bins, followed by merging fragmented bins and recruiting unbinned contigs.

## Install and Setup
### conda
We recommend using conda to run HiFine.

After installing Anaconda (or miniconda), Users can clone the repository with git:
```
git clone https://github.com/dyxstat/HiFine.git
```

Once complete, you can enter the repository folder and then create a HiFine environment using conda:
```
# enter the HiFine folder
cd HiFine
# Construct environment
conda env create -f HiFine_linux_env.yaml 
or
conda env create -f HiFine_osx_env.yaml
# Enter the environment
conda activate HiFine_env
```
