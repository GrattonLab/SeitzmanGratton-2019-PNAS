function findVariantSubgroups(varCorrLoc)
%findVariantSubgroups(varCorrLoc)
%
% This function clusters individuals into subgroups on the basis of the
% correlation coefficient between the individual's network variants and all
% template networks. Must be run after templateMatchingVariants.m.
%
% INPUT
%
% varCorrLoc: a path to the correlation coefficient files generated from
% template_matching_variants.m (the files that end in corrCoeffFor*).
%
% OPTIONAL INPUT
%
% outputdir: the directory to which the output files will be written
% 
% OUTPUTS
% two figures: the first is the sorted subject-subject correlation matrix,
% and the second is the average (across all subjects) of the mean 
% similarity to the template networks (as in Figure 5 from Seitzman et al.,
% 2019).
%
% "Where there's a will there's a kluge."
% -BAS 10/11/2019

%%%%% CHANGE THIS PATH TO THE LOCATION WHERE YOU STORED THE TEMPLATE NAMES %%%%%

% templateLoc = '/your/path/here/';
templateLoc = pwd;
load([templateLoc '/networkTemplateNames.mat']);
count = 0;
for i=1:length(IDNames)
    if strcmp(IDNames{i},'skip')
        continue;
    else
        count = count + 1;
        tempIDNames{count} = IDNames{i};
    end
end
IDNames = tempIDNames;
%%%%% CHANGE THIS PATH TO THE LOCATION WHERE YOU STORED THE TEMPLATE NAMES %%%%%


% Set variables
if ~exist('outputdir')
    outputdir = pwd;
end


% Determine number of subjects
% fileNameStart = [varCorrLoc '/variantTemplatematch_allSubjects_corrCoeffFor'];
fileNameStart = [varCorrLoc '/Templatematch_spCorr_bysubject_corrCoeffFor'];
tempFile = ft_read_cifti_mod([fileNameStart IDNames{1} '.dtseries.nii']);
numSubs = size(tempFile.data,2);
clear tempFile
meanCorr = zeros(length(IDNames),numSubs); 


% Compute variant-size-weighted-average of all variants' correlation 
% coefficient to each template network (average across all variants within
% each subject)
for i=1:length(IDNames) 
    dataTemp = ft_read_cifti_mod([fileNameStart IDNames{i} '.dtseries.nii']);
    dataTemp = dataTemp.data;
    dataTemp(dataTemp==0) = NaN;
    meanCorr(i,:) = nanmean(dataTemp); 
end
clear dataTemp


% Correlate this measure across all individuals
subBySubCorrMat = corr(meanCorr); % should be NxN where N = numSubs


% Cluster across individuals (use your preferred clustering algorithm here)
subgroupIDs = cluster(linkage(subBySubCorrMat,'weighted'),'maxclust',2);
subGroups = unique(subgroupIDs);
for i=1:length(subGroups)
    groupInds{i} = find(subgroupIDs==subGroups(i));
end
[~,sortedInds] = sort(subgroupIDs);


% Plot the results
% sorted matrix
figure;
imagesc(subBySubCorrMat(sortedInds,sortedInds),[-1 1]); colorbar; colormap(redblue)
set(gca,'FontName','Arial','FontSize',16)

% average similarity
netOrder = [8 9 10 2 4 7 3 5 1 12 6 13 14 11]; % from Seitzman et al., 2019 (rearrange as desired)
figure;
hold on
for i=1:length(subGroups)
    groupMeanVals = meanCorr(netOrder,groupInds{i})';
    subgroupColor = rand(1,3);
    errorbar(1:length(IDNames),mean(groupMeanVals),std(groupMeanVals)./sqrt(length(groupInds{i})),'-','Color',subgroupColor,'LineWidth',2)
end
ylabel('mean similarity to network template','FontName','Arial','FontSize',16)
set(gca,'FontName','Arial','FontSize',16,'xtick',1:length(IDNames),'xticklabel',IDNames(netOrder))


end

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 

end