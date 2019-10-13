function createSpatialCorrMap(groupDconnLoc,subDconnLoc,outputdir)
%createSpatialCorrMap(dconnLoc,subDataLoc,[outputdir])
%
% This function creates a spatial correlation map by comparing cortical
% BOLD data from a single individual to a group-average. The assumed file
% format is CIFTI, in order to work with Connectome Workbench (see 
% https://www.nitrc.org/projects/cifti/ for more details). This code 
% requires the GIFTI and CIFTI_Resources packages to be added to the user's
% path (which released in addition to this one).
%
% INPUTS
% dconnLoc: a path to the group-average correlation matrix (dconn is CIFTI
% for correlation matrix) and the file name
%
% subDataLoc: a path to the single individual's correlation matrix (in the 
% CIFTI format) and the file name
%
% OPTIONAL INPUT
% outputdir: the directory to which the output file will be written
% 
% OUTPUT
% a single CIFTI file that contains the spatial correlation map 
%
% "Where there's a will there's a kluge."
% -BAS 10/11/2019


%%%%% CHANGE THIS PATH TO THE LOCATION WHERE YOU STORED THE TEMPLATE %%%%%
%%%%% NOTE: This template will only work with the 32k-fsLR surfaces. %%%%%

% templateLoc = '/your/path/here/';
templateLoc = pwd;
template = ft_read_cifti_mod([templateLoc '/templateSpatialCorrMap.dtseries.nii']);
template.data = [];

%%%%% CHANGE THIS PATH TO THE LOCATION WHERE YOU STORED THE TEMPLATE %%%%%


% Set variables
if ~exist('outputdir')
    outputdir = pwd;
end


% Read in group-average matrix
tempCifti = ft_read_cifti_mod(groupDconnLoc);
cortexInds = 1:sum(tempCifti.brainstructure==1 | tempCifti.brainstructure==2);
groupMat = single(FisherTransform(tempCifti.data(cortexInds,cortexInds)));
clear tempCifti


% Read in single subject matrix
tempCifti = ft_read_cifti_mod(subDconnLoc);
cortexInds = 1:sum(tempCifti.brainstructure==1 | tempCifti.brainstructure==2);
subMat = single(FisherTransform(tempCifti.data(cortexInds,cortexInds)));
clear tempCifti


% Compare single subject to group-average at each cortical location
for i=1:length(cortexInds)
    template.data(i,1) = paircorr_mod(groupMat(:,i),subMat(:,i));
end
clear groupMat subMat


% Write out the spatial correlation map
ft_write_cifti_mod([outputdir '/spatialCorrMap.dtseries.nii'],template)

end

