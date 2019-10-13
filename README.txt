The included MATLAB scripts were used in the study: Seitzman et al., 2019.
Trait-like variants in human functional brain networks. PNAS. These scripts
will allow users to create network variants, assign them to network 
templates, and find subgroups/sub-types of individuals on the basis of
network variants. The order in which the scripts should be run are:

1. createSpatialCorrMap.m
2. binarizeAndIDvariants.m
3. templateMatchingVariants.m
(optional 3.5. reassignVariants.m)
4. clusterIndividuals.m
5. findVariantSubgroups.m

All of these scripts require supporting scripts found in the Resources 
folder to be a part of your MATLAB path. One option is to run the following
command in the MATLAB command window before running the first script:

addpath(genpath('/your/path/to/the/Resources/folder/Resources/'))

Some of the scripts require supporting .mat or .dtseries.nii (CIFTI) files,
which are also included in this release. 

Please cite your usage of these scripts (Seitzman et al., 2019. Trait-like 
variants in human functional brain networks. PNAS). For any questions, 
problems, or bugs found, please email seitzman@wustl.edu and/or 
cgratton@northwestern.edu. Thank you!