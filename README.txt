The included MATLAB scripts were used in the study: Seitzman et al., 2019.
Trait-like variants in human functional brain networks. PNAS. These scripts
will allow users to create network variants, assign them to network 
templates, and find subgroups/sub-types of individuals on the basis of
network variants. The order in which the scripts should be run are:

1. createSpatialCorrMap.m
2. binarizeAndIDvariants.m
3. templateMatchingVariants.m
(optional 3.5. reassignTemplates.m)
4. findVariantSubgroups.m

All of these scripts require supporting scripts for reading and writing 
CIFTI files (see more at https://github.com/fieldtrip/fieldtrip/). Modified
version of these scripts are provided in the Resources folder. This folder 
also includes two scripts for computing pairwise correlations in parallel 
and Fisher-Z transforming correlation coefficients. One option is to run 
the following command in the MATLAB command window before running the first
script: addpath(genpath('/your/path/to/the/Resources/folder/Resources/'))

Additionally, supporting .mat or .dtseries.nii (CIFTI) files are included 
in this release. Please cite our paper when using any of these scripts  
(Seitzman et al., 2019. Trait-like variants in human functional brain 
networks. PNAS). For any questions, problems, or bugs found, please email 
seitzman@wustl.edu and/or cgratton@northwestern.edu. Thank you!