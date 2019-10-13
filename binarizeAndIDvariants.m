function binarizeAndIDvariants(spatialCorrMap,snrMask,corrThresh,sizeThresh,outputdir)
%binarizeAndIDvariants(spatialCorrMap,[snrMask],[corrThresh],[sizeThresh],[outputdir])
%
% This function creates network variants from a spatial correlation map and
% gives each contiguous variant (of sufficient size) a unique ID. Must be
% run after createSpatialCorrMap.m. This code requires a neighbors file,
% which indicates the numbers of all neighboring vertices for each cortical
% vertex (again, CIFTI format is assumed).
% 
% INPUT
% spatialCorrMap: a path to the spatial correlation map (generated in the
% previous step) and the file name
%
% OPTIONAL INPUTS
% snrMask: the path and file name of a CIFTI containing a map of regions 
% to be excluded on the basis of some signal-to-noise measure. The map is 
% assumed to be binary (0 or 1). If a variant touches any part of the map
% where there is 1, it will be excluded from all further analyses.
%
% variantThresh: the threshold below which vertices will be considered
% network variants (all verts < x set to 1, all others set to 0). If a
% value is not specified, the threshold will be set to 0.1 (lowest decile).
% 
% sizeThresh: the threshold below which network variants will be excluded
% (variants with N contiguous verts < x will be set to 0). If a value is
% not specified, the threshold will be set to 30.
%
% outputdir: the directory to which the output files will be written
% 
% OUTPUTS
% two CIFTI files, the first contains binary network variants and the other
% contains unique ID labels for each network variant
%
% "Where there's a will there's a kluge."
% -BAS 10/11/2019


%%%%% CHANGE THIS PATH TO THE LOCATION WHERE YOU STORED THE NEIGHBORS %%%%%
%%%%% NOTE: This neighbors file is specific to the 32k-fsLR surfaces. %%%%%
% This file is a 59412x7 matrix of neighboring surface vertices. For each
% surface vertex (first column), it's hexagonal neighboring vertices are
% listed (columns 2-7). If no neighbor exists (e.g., medial wall), the
% column will display NaN. WARNING: If you are not using the Conte69 32k-
% fsLR atlas surfaces, a new neighbors file must be made. 

% neighLoc = '/your/path/here/';
neighLoc = pwd;
neigh = load([neighLoc '/Cifti_surf_neighbors_LR_normalwall.mat']);
neigh = neigh.neighbors;

%%%%% CHANGE THIS PATH TO THE LOCATION WHERE YOU STORED THE NEIGHBORS %%%%%


% Set variables
template = ft_read_cifti_mod(spatialCorrMap);
corrData = template.data;
cortexInds = 1:sum(template.brainstructure==1 | template.brainstructure==2);


if ~exist('snrMask')
    mask = zeros(length(cortexInds),1);
    mask = logical(mask==1);
else
    mask = ft_read_cifti_mod(snrMask);
    mask = logical(mask.data(cortexInds)==1);
end


if ~exist('corrThresh')
    cThresh = 0.1;
else
    cThresh = corrThresh;
end
variantThresh = 100*cThresh;

if ~exist('sizeThresh')
    sThresh = 30;
else
    sThresh = sizeThresh;
end

if ~exist('outputdir')
    outputdir = pwd;
end


% Binarize vertices (1 = low spatial correlation value = potential variant)
sorted = sort(corrData);
cThresh = sorted(floor(cThresh*length(sorted)));
clear sorted
binData = logical(corrData<cThresh);
template.data = binData;


% Write out binary variants (no size threshold applied)
ft_write_cifti_mod([outputdir '/binarySpatialCorrMap_thresh' num2str(variantThresh) '%.dtseries.nii'],template)


% Give each variant a unique ID.
netVars = zeros(length(cortexInds),1);
count=1;
for i=cortexInds
    
    % If the vertex is a potential variant,
    if binData(i)==1
        
        % And if the vertex is already part of a variant, use existing ID
        if netVars(i)>0
            id = netVars(i);

        % Otherwise, use a new ID
        else
            id = count;
            count = count+1;
        end
        
        % Label the vertex
        netVars(i) = id;
        
        % Check that vertex's neighbors
        neighVerts = neigh(i,2:7);
        neighVerts(isnan(neighVerts)) = []; % Ignore NaNs (medial wall vertices)
        
        % All nonzero neighbors are given the same ID
        for j=1:length(neighVerts)
            % Do not overwrite previously assigned variants
            if netVars(neighVerts(j))>0
                netVars(netVars==netVars(neighVerts(j)))=id;
            elseif binData(neighVerts(j))==1 && netVars(neighVerts(j))==0
                netVars(neighVerts(j))=id;
            end
        end
    end
    
end


% Elimiate variants that are too small
ids = unique(netVars); 
if ids(1)==0; ids(1) = []; end; % Zero is not a variant label
for i=1:length(ids)
    inds = find(netVars==ids(i));
    if length(inds)<sThresh
        netVars(inds)=0;
    end
end


% Eliminate variants that touch the SNR mask
ids = unique(netVars); 
if ids(1)==0; ids(1) = []; end; % Zero is not a variant label
for i=1:length(ids)
    inds = find(netVars==ids(i));
    if sum(mask(inds))>0
        netVars(inds)=0;
    end
end


% Normalize IDs (to be sequential)
ids = unique(netVars);
if ids(1)==0; ids(1) = []; end; % Zero is not a variant label
for i=1:length(ids)
    netVars(logical(netVars==ids(i)))=i;
end


% Write out variants with uniqueIDs
template.data=netVars;
ft_write_cifti_mod([outputdir '/networkVariants_thresh' num2str(variantThresh) '%_sizeThresh' num2str(sThresh) '.dtseries.nii'],template)

end