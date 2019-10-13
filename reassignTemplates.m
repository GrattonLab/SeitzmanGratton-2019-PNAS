function reassignTemplates(varList,varCorrLoc,varAssignLoc,groupAssign,outputdir)
%reassignTemplates(varList,varCorrLoc,varAssign,[groupAssign],[outputdir])
%
% This function re-assigns network labels to variants after considering the
% strength of the winning template network and the spatial overlap between 
% the variant (with network assignment) and the template network. Variants
% that have a winning template network correlation coefficient < 0.3 are
% assigned to "unknown" network. Variants with a spatial overlap of >50%
% between the variant and the same group-average template network are
% removed (e.g., if a variant is assigned to Default Mode and that variant
% spatially overlaps with the template Default Mode Network by >50%, it is
% removed). The function is designed to run as a single batch for all 
% subjects. Must be run after templateMatchingVariants.m.
%
% INPUTS
% varList: a path to the location of the network variants for each subject
% (and the file name). This is the file that contains a unique ID for each
% network variant.
%
% varCorrLoc: a path to the correlation coefficient files generated from
% template_matching_variants.m (the files that end in corrCoeffFor*).
%
% varAssignLoc: a path to the networkIDs file generated from
% template_matching_variants.m.
%
% OPTIONAL INPUTS
% groupAssign: a path to the group-average networks used as the templates
% (and the file name)
%
% outputdir: the directory to which the output files will be written
% 
% OUTPUTS
% a single CIFTI that contains the network label of each variant for each
% subject. Each "time point" (column) in the CIFTI corresponds to one
% subject. The order of subjects will match the order of the input files
% (which are specified in the params_file).
%
% a .mat file containing the number of variants for each subject before and
% after running this code
%
% "Where there's a will there's a kluge."
% -BAS 10/11/2019


%%%%% CHANGE THIS PATH TO THE LOCATION WHERE YOU STORED THE TEMPLATE NAMES %%%%%

% templateLoc = '/your/path/here/';
templateLoc = pwd;
load([templateLoc '/networkTemplateNames.mat']);

%%%%% CHANGE THIS PATH TO THE LOCATION WHERE YOU STORED THE TEMPLATE NAMES %%%%%


% Set variables
if ~exist('outputdir')
    outputdir = pwd;
end

if ~exist('groupAssign')  
    %%%%% CHANGE THIS PATH TO THE LOCATION WHERE YOU STORED THE GROUP-AVERAGE NETWORKS %%%%%
    
    %groupLoc = '/your/path/here/';
    groupLoc = pwd;
    group = ft_read_cifti_mod([groupLoc '/WashU120_groupNetworks.dtseries.nii']);
    
    %%%%% CHANGE THIS PATH TO THE LOCATION WHERE YOU STORED THE GROUP-AVERAGE NETWORKS %%%%%  
else
    group = ft_read_cifti_mod(groupAssign);
end


% Read in location of variants
varLoc = textread(varList,'%s');
numVars = zeros(size(varLoc,1),2);


% Load variants' network assignments
varAssn = ft_read_cifti_mod([varAssignLoc '/variantTemplatematch_allSubjects_networkIDs.dtseries.nii']);
varNet = varAssn.data;
cortexInds = 1:sum(varAssn.brainstructure==1 | varAssn.brainstructure==2);
group = group.data(cortexInds,1);


% Load variants' correlation coefficients to each template network
numNetsCount = 0;
for i=1:length(IDNames)
    if strcmp(IDNames{i},'skip')
        varCorrsAll{i} = 0;
    else
        numNetsCount = numNetsCount + 1;
        varCorrsAll{i} = ft_read_cifti_mod([varCorrLoc '/variantTemplatematch_allSubjects_corrCoeffFor' IDNames{i} '.dtseries.nii']);
    end
end


% Loop through subjects
for i=1:length(varLoc)
    
    
    % Load variants (uniqueIDs)
    disp(['Re-assigning variants for subject ' num2str(i)])
    vars = ft_read_cifti_mod(varLoc{i});
    varsData = vars.data; 
    varIDs = unique(varsData); 
    if varIDs(1)==0; varIDs(1)=[]; end;
    
    
    % Store original number of variants
    numVars(i,1) = length(varIDs);
    
    
    % Look for variants with poor match to the winning template network
    matches = zeros(length(varIDs),1);
    reassign = zeros(length(varIDs),1);
    for j=1:length(varIDs)
        % Find the variant's network assignment
        inds = find(varsData==varIDs(j));
        indivNet = varNet(inds(1),i);
        
        
        % Extract the correlation coefficient (i.e., the match to that
        % template network)
        for qq=1:numNetsCount
            if indivNet==qq
                varCorrs = varCorrsAll{qq}.data;
                varCorrs = varCorrs(inds(1),i);
                break
            end
        end
        
        
        % If it's a poor match (corr coeff<0.3), flag the variant for 
        % reassignment
        if varCorrs<0.3
            reassign(j)=1;
            continue
        end
        
        
        % Check the group-average networks at that location
        groupNet = group(inds,1);
        numNet = unique(groupNet);
        
        
        % If the variant's network matches the group-average network, and
        % the spatial overlap between them is >50%, flag that variant.
        if length(numNet)>1
            if ismember(indivNet,groupNet)
                ind = zeros(length(numNet),1);
                for k=1:length(numNet)
                    ind(k) = length(find(groupNet==numNet(k)));
                end
                ind = ind./length(groupNet);
                [val,loc] = max(ind);
                if val>0.5 && numNet(loc)==indivNet
                    matches(j)=1;
                end
            else
                continue
            end
        else
            if indivNet==groupNet(1)
                matches(j)=1;
            else
                continue
            end
        end
    end
    
    
    % Re-assign or remove flagged variants
    for j=1:length(varIDs)
        if reassign(j)==1
            varNet(logical(varsData==varIDs(j)),i)=6; % "unknown" network
        elseif matches(j)==1
            varNet(logical(varsData==varIDs(j)),i)=0; % remove
        end
    end
    
    
    % Store the new number of variants (after removal of matches)
    numVars(i,2) = length(varIDs)-nnz(matches);
end


% Write out the results
save('numVariants_beforeAfterReassign.mat','numVars')
varAssn.data = varNet;
ft_write_cifti_mod([outputdir '/variantTemplatematch_allSubjects_networkIDs_reassigned.dtseries.nii'],varAssn)

end