# BPPmulti
Perl script for running multiple iBPP and BPP instances

## Description
Automatically run iBPPv2.1.2 (Solís-Lemus et al. 2015; https://github.com/cecileane/iBPP/) and BPPv3 (Yang & Rannala 2014; http://abacus.gene.ucl.ac.uk/software/) multiple times
- with pairwise combinations of X tau and Y theta priors
- on multiple data sets (prior, molecular, trait, integrative)
- on multiple guide tree topologies
- in parallel

## Usage
Directly edit the 'User Settings' section in the perl script (also see comments within the script):
- Fill in the paths to your working directory and to your files. The files should have the format used by (i)BPP.
- Set your options. They are used to make the BPP config-files.
- Enter your guide tree topologies.
- Switch your desired tasks on/off.
- Then save the file and execute it with perl. Make the file executable on Unix systems. Subfolders will be created in your working directory for each analysis. The results will be collected to 'final_output.txt' in your working directory.


Windows: >perl BPPmulti.pl

Unix:    >./BPPmulti.pl

## References
Solís-Lemus, C., L. L. Knowles, and C. Ané. 2015. Bayesian species delimitation combining multiple genes and traits in a unified framework. Evolution (N. Y). 69:492–507.

Yang, Z., and B. Rannala. 2014. Unguided Species Delimitation Using DNA Sequence Data from Multiple Loci. Mol. Biol. Evol. 31:3125–3135.
