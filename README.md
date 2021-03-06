# HiFine: integrating Hi-c-based and shotgun-based methods to reFine binning of metagenomic contigs

## Introduction
HiFine is a pipeline to refine the binning results of metagenomic contigs by integrating both Hi-C-based and shotgun-based binning tools. HiFine designs a strategy of fragmentation for the original bin sets derived from the Hi-C-based and shotgun-based binning methods, which considerably increases the purity of initial bins, followed by merging fragmented bins and recruiting unbinned contigs.

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

## Preparation
### Deriving the initial sets from different types of binning methods
Two initial binning sets are needed, corresponding to two folders. One folder contains bins generated by a Hi-C-based binning method and the second folder includes bins constructed by a shotgun-based binning pipeline.

### Generating metagenomic Hi-C contact maps normalized by HiCzin
* If you employs [HiCBin](https://github.com/dyxstat/HiCBin) to bin the contigs, you can directly get access to normalized Hi-C contact maps through the file *‘HiCzin\_normalized\_contact.gz’* in the output directory of HiCBin.
* If you use other Hi-C-based binning methods such as ProxiMeta or [bin3C](https://github.com/cerebis/bin3C), 
you need to run the [**HiCzin**](https://github.com/dyxstat/HiCzin) software (https://github.com/dyxstat/HiCzin) to generate normalized Hi-C contact maps. The normalized Hi-C contacts and other needed contig information are stored in the
file *'HiCzin\_normalized\_contact.gz'* from the output directory of HiCzin software.

**Then you can move the file *'HiCzin\_normalized\_contact.gz'* to the directory of HiFine folder.**


## Usage
### Implement the HiFine pipeline 
```
python ./hifine.py refine [Parameters] HiC_folder Shotgun_folder HiC_contacts Contig_file Output_directory
```
#### Parameters
```
--min-binsize: Minimum bin size of shared bins [default: 500000]
--min-complete: Minimum bin size of relatively complete bins [default: 500000]
--min-frac: Fraction to determine the relatively complete bins [default: 0.8]
--alpha: hyperparameter in step2 [default: 0.6]
--beta: hyperparameter in step3 [default: 0.3]
--bin3C: Whether the Hi-C-based binning method is bin3C [default: False]
--cover: Cover existing files
-v: Verbose output
```

#### Input File

* *HiC_folder*: a folder path containing bins generated by a Hi-C-based binning method
* *Shotgun_folder*: a folder path including bins constructed by a shotgun-based binning method
* *HiC_contacts*: Hi-C contact maps normalized by HiCzin (i.e., file 'HiCzin_normalized_contact.gz')
* *Contig_file*: a fasta file of the assembled contig (e.g. final.contigs.fa)


#### Example
```
python ./hifine.py refine -v hicbin_bin_folder metabat2_bin_folder HiCzin_normalized_contact.gz final.contigs.fa out
```
Specically, if you employs bin3C as the Hi-C-based binning method, use
```
python ./hifine.py refine --bin3C -v bin3C_bin_folder metabat2_bin_folder HiCzin_normalized_contact.gz final.contigs.fa out
```
The results of the pipeline action are all in the 'out' directory.

#### Output File
* *HIFINE_BIN:* Final refined genomic bins by HiFine
* *hifine_cluster.txt:* Specific group labels of binned contigs
* *hifine.log:* log file of the specific implementation information of HiCBin


## Contacts and bug reports
If you have any questions or suggestions, welcome to contact Yuxuan Du (yuxuandu@usc.edu).


## Copyright and License Information
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.






