function templateMatchingVariants(paramsFile,outputdir)
%templateMatchingVariants('paramsFile.m',[outputdir])
%
% This function assigns network labels to variants by assessing the 
% similarity of each variant's connectivity map to a set of template
% network connectivity maps. Each variant is assigned the identity of the 
% network with the most similar template network connectivity map. The
% function is designed to run as a single batch for all subjects. Must be
% run after binarizeAndIDvariants.m. Modified from Evan Gordon's original
% template matching code (see Gordon et al., 2017. Individual-specific 
% features of brain systems identified with resting state functional 
% correlations. NeuroImage).
%
% INPUT
% params_file: a parameters file (a .m file) which will be executed to load needed 
% parameters, including:
% - a datalist and tmasklist for the subjects to be run
% - a set of template network connectivity maps
%
% OPTIONAL INPUT
% outputdir: the directory to which the output files will be written
% 
% OUTPUTS
% a single CIFTI that contains the network label of each variant for each
% subject. Each "time point" (column) in the CIFTI corresponds to one
% subject. The order of subjects will match the order of the input files
% (which are specified in the params_file).
%
% a single CIFTI that contains the ratio of the top two template matches
% for each variant. Each "time point" (column) in the CIFTI corresponds to 
% one subject. The order of subjects will match the order of the input 
% files (which are specified in the params_file).
%
% a number of CIFTIs that contain the correlation coefficient between each
% variant and the specified template network. There will be one file per
% template network. Each "time point" (column) in the CIFTI corresponds to 
% one subject. The order of subjects will match the order of the input 
% files (which are specified in the params_file).
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


% Load variables specified in the params_file
[paramspath,paramsname,paramsextension] = fileparts(paramsFile);
origpath = pwd;
if ~isempty(paramspath)
    cd(paramspath)
end
params = feval(paramsname);
varnames = fieldnames(params);
for i = 1:length(varnames)
    evalc([varnames{i} ' = params.' varnames{i}]);
end
clear varnames params
cd(origpath)


% Load subject lists, data files, and network templates
[subjects, ciftifiles] = textread(surfdatafile,'%s %s');
[subjects, tmasks] = textread(tmaskfile,'%s %s');
variants = textread(variantsfile,'%s');
load(templatesfile)
ThreshTemplates = templates;
clear templates


% Loop through subjects
count=0;
prevstring = [];
for s = 1:length(subjects)
    count=count+1;
    subject = subjects{s};
    
    % Load timecourse data and variant information
    string = ['Subject ' subject ': calculating variant correlation maps'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    variant_file = variants{count};
    cifti_file = ciftifiles{s};
    tmask = dlmread(tmasks{s});
    subdata = ft_read_cifti_mod(cifti_file);
    subdata = subdata.data;
    
    netVars = ft_read_cifti_mod(variant_file);
    uniqueVars = unique(netVars.data); 
    if uniqueVars(1)==0; uniqueVars(1)=[]; end;
    

    % Set up needed variables if it's the first iteration
    if s==1;
        out_template = netVars; out_template.data = [];
        ncortverts = sum(netVars.brainstructure==1 | netVars.brainstructure==2);
        networkIDs = zeros(ncortverts,length(subjects));
        networkRatio = zeros(ncortverts,length(subjects));
        networkPercent = zeros(ncortverts,size(ThreshTemplates,2),length(subjects));
    end
    
    subdata = subdata(1:ncortverts,logical(tmask));

    
    % Loop through each variant and calculate its connectivity map
    for jj=1:length(uniqueVars)
        varverts = logical(netVars.data==uniqueVars(jj));
        correlmaps = paircorr_mod(mean(subdata(varverts,:))',subdata'); 
        correlmaps(isnan(correlmaps)) = 0;
        correlmaps_var(:,jj) = correlmaps';
        clear correlmaps
    end
    disp('Done')
    clear subdata
    
    
    % Match each variant's connectivity map to each template connectivity map
    disp('Matching to templates and computing goodness of fit')
    corr_coeff = zeros(length(uniqueVars),size(ThreshTemplates,2));
    for var=1:length(uniqueVars)
        disp(['Subject ' subject ', variant #' num2str(var) ' out of ' num2str(length(uniqueVars))]);
        for templatenum = 1:size(ThreshTemplates,2);
            corr_coeff(var,templatenum) = paircorr_mod(correlmaps_var(:,var),ThreshTemplates(:,templatenum));
        end
    end
    
    
    % Determine the 1st and 2nd place template network (max and next max match)
    corr_coeff(isnan(corr_coeff)) = 0;
    [~,maxi] = max(corr_coeff,[],2);
    tempCorrCoeff = corr_coeff;
    for qq=1:length(uniqueVars)
        tempCorrCoeff(qq,maxi(qq))=0;
    end
    [~,nextMax] = max(tempCorrCoeff,[],2);
    clear tempCorrCoeff
    
    
    % Save the correlation between each variant and each template network
    % (between their connectivity maps)
    for var=1:length(uniqueVars)
        varverts = logical(netVars.data==uniqueVars(var));
        corrVals = corr_coeff(var,:);
        networkIDs(varverts,count) = IDs(maxi(var));
        networkRatio(varverts,count) = corr_coeff(var,maxi(var))./corr_coeff(var,nextMax(var));
        for kk=1:size(ThreshTemplates,2)
            networkPercent(varverts,kk,count)=corrVals(kk);
        end
    end
    clear maxi corr_coeff nextMax corrVals
end


% Write out the results
out_template.data = zeros(ncortverts,count);

out_template.data(1:size(networkIDs,1),:) = networkIDs;
ft_write_cifti_mod([outputdir '/variantTemplatematch_allSubjects_networkIDs'],out_template);
clear networkIDs

out_template.data(1:size(networkRatio,1),:) = networkRatio;
ft_write_cifti_mod([outputdir '/variantTemplatematch_allSubjects_ratioOfTopTwoTemplates'],out_template);
clear networkRatio

for kk = 1:size(ThreshTemplates,2)
    out_template.data(1:size(networkPercent,1),:) = squeeze(networkPercent(:,kk,:));
    ft_write_cifti_mod([outputdir '/variantTemplatematch_allSubjects_corrCoeffFor' IDNames{kk}],out_template);
end

end

function [pathstr, name, ext] = fileparts(file)
%FILEPARTS Filename parts.
%   [PATHSTR,NAME,EXT] = FILEPARTS(FILE) returns the path, file name, and
%   file name extension for the specified FILE. The FILE input is a string
%   containing the name of a file or folder, and can include a path and
%   file name extension. The function interprets all characters following
%   the right-most path delimiter as a file name plus extension.
%
%   If the FILE input consists of a folder name only, be sure that the
%   right-most character is a path delimiter (/ or \). Othewise, FILEPARTS
%   parses the trailing portion of FILE as the name of a file and returns
%   it in NAME instead of in PATHSTR.
%
%   FILEPARTS only parses file names. It does not verify that the file or
%   folder exists. You can reconstruct the file from the parts using
%      fullfile(pathstr,[name ext])
%
%   FILEPARTS is platform dependent.
%
%   On Microsoft Windows systems, you can use either forward (/) or back
%   (\) slashes as path delimiters, even within the same string. On Unix
%   and Macintosh systems, use only / as a delimiter.
%
%   See also FULLFILE, PATHSEP, FILESEP.

%   Copyright 1984-2012 The MathWorks, Inc.
%   $Revision: 1.18.4.18 $ $Date: 2012/04/14 04:15:41 $

pathstr = '';
name = '';
ext = '';

if ~ischar(file)
    error(message('MATLAB:fileparts:MustBeChar'));
elseif isempty(file) % isrow('') returns false, do this check first
    return;
elseif ~isrow(file)
    error(message('MATLAB:fileparts:MustBeChar'));
end

if ispc
    ind = find(file == '/'|file == '\', 1, 'last');
    if isempty(ind)
        ind = find(file == ':', 1, 'last');
        if ~isempty(ind)       
            pathstr = file(1:ind);
        end
    else
        if ind == 2 && (file(1) == '\' || file(1) == '/')
            %special case for UNC server
            pathstr =  file;
            ind = length(file);
        else 
            pathstr = file(1:ind-1);
        end
    end
    if isempty(ind)       
        name = file;
    else
        if ~isempty(pathstr) && pathstr(end)==':' && ...
                (length(pathstr)>2 || (length(file) >=3 && file(3) == '\'))
                %don't append to D: like which is volume path on windows
            pathstr = [pathstr '\'];
        elseif isempty(deblank(pathstr))
            pathstr = '\';
        end
        name = file(ind+1:end);
    end
else    % UNIX
    ind = find(file == '/', 1, 'last');
    if isempty(ind)
        name = file;
    else
        pathstr = file(1:ind-1); 

        % Do not forget to add filesep when in the root filesystem
        if isempty(deblank(pathstr))
            pathstr = '/';
        end
        name = file(ind+1:end);
    end
end

if isempty(name)
    return;
end

% Look for EXTENSION part
ind = find(name == '.', 1, 'last');

if isempty(ind)
    return;
else
    ext = name(ind:end);
    name(ind:end) = [];
end

end

