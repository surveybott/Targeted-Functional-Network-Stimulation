function pfm_ninet(sub,varargin)
% inputs
tmpDir = getenv('TMPDIR');
p = inputParser;
p.addRequired('sub');
p.addParameter('ses',[]);
p.addParameter('derivativeDir','/scratch/st-fidelvil-1/LeRNIT/bids/derivatives');
p.addParameter('workingDir','/scratch/st-fidelvil-1/LeRNIT/working/pfm/');
p.addParameter('softwareDir','/arc/project/st-fidelvil-1/software');
p.addParameter('funcSuffix','desc-optcomDenoised_bold*.nii')
p.addParameter('distMatrix','/arc/project/st-fidelvil-1/software/Targeted-Functional-Network-Stimulation/PFM/DistanceMatrix.mat');
p.addParameter('spatialThresh',[20 20],@(x) isnumeric(x) && numel(x)==2); % [CortexThreshold, SubcortexThreshold] for spatial_filtering.m
p.addParameter('clusterDir',tmpDir);
p.parse(sub,varargin{:});
inputs = p.Results;
if ~isempty(inputs.ses) 
    if ischar(inputs.ses)
        inputs.ses = {inputs.ses};
    end
    inputs.ses = regexprep(inputs.ses,'ses-','');
end
% add code
code = {'Targeted-Functional-Network-Stimulation', 'MSCcodebase'};
for i=1:numel(code)
    addpath(genpath(fullfile(inputs.softwareDir,code{i})));
end
% setup parallel
if ~exist('cores','var')
    cores = feature('numcores') - 1;
end
pc = parcluster;
if ~isempty(inputs.clusterDir)
    pc.JobStorageLocation = inputs.clusterDir;
end
warning('off','stats:regress:RankDefDesignMat');
% setup data/paths/ get sessions
sessions = dir(fullfile(inputs.derivativeDir,['sub-' sub],'ses-*'));
ses = regexprep({sessions.name},'ses-','');
if ~isempty(inputs.ses)
    ses = ses(ismember(ses,inputs.ses));
    if isempty(ses)
        error('No sessioins found');
    end
end
for i=1:numel(ses)
    tmpDir = fullfile(inputs.workingDir,sprintf('sub-%s_ses-%s',sub,ses{i}));
    outDir = fullfile(inputs.derivativeDir,['sub-' sub],['ses-' ses{i}],'func');
    if ~isfolder(tmpDir)
        mkdir(tmpDir);
    end
    fprintf('sub-%s_ses-%s\n',sub,ses{i});

    % ciftis
    tmp = dir(fullfile(outDir,['*' inputs.funcSuffix]));
    cifti = arrayfun(@(x) fullfile(x.folder,x.name),tmp,'UniformOutput',0);
    tmp = dir(fullfile(inputs.softwareDir,'/templateflow/tpl-fsLR/tpl-fsLR_den-32k_*_midthickness.surf.gii'));
    surf = arrayfun(@(x) fullfile(x.folder,x.name),tmp,'UniformOutput',0);

    % create/load distance matrix
    distMat = inputs.distMatrix;
    if ~exist(distMat,'file')
        make_distance_matrix(cifti{1},surf,tmpDir,{pc cores});
        distMat = fullfile(tmpDir,'DistanceMatrix.mat');
    end

    % load cifti, do global signal regression, concatenate
    c_concat = [];
    for j=1:numel(cifti)
        c = ft_read_cifti_mod(cifti{j});
        % GSR
        cortex_idx = find(c.brainstructure == 1 | c.brainstructure == 2);
        gs = nanmean(c.data(cortex_idx,:),1);
        for v=cortex_idx'
            [~,~,c.data(v,:)] = regress(c.data(v,:)',[gs' ones(size(c.data,2),1)]);
        end
        c.data = c.data - mean(c.data,2); %demean
        c_concat = [c_concat c.data];
    end
    c.data = c_concat;
    c.hdr.dim(6) = size(c.data,2);
    clear c_concat

    % subcortical regression
    c = regress_cortical_signals(c, distMat, 20);
    cifti_reg = fullfile(tmpDir,sprintf('sub-%s_ses-%s_task-restME_space-fsLR_den-91k_desc-optcomDenoisedRegressed_bold.dtseries.nii',sub,ses{i}));
    ft_write_cifti_mod(cifti_reg,c);


    % smooth with geodesic (for surface data) and Euclidean (for volumetric data) Gaussian kernels;
    kernel.cort = 2.55;
    kernel.subcort = 2.55;
    cifti_smooth = regexprep(cifti_reg,'DenoisedRegressed','DenoisedRegressedSmoothed');
    system(['wb_command -cifti-smoothing ' cifti_reg ' ' num2str(kernel.cort) ' ' num2str(kernel.subcort) ' COLUMN ' cifti_smooth ' -left-surface ' surf{1} ' -right-surface ' surf{2} ' -merged-volume']);

    % load concat/preprocessed cifti and FD
    c = ft_read_cifti_mod(cifti_smooth);
    tsv = dir(sprintf('/scratch/st-fidelvil-1/LeRNIT/bids/derivatives/sub-%s/ses-%s/func/*task-restME_*desc-confounds_timeseries.tsv',sub,ses{i}));
    fd = [];
    for j=1:numel(tsv)
        t = tdfread(fullfile(tsv(j).folder,tsv(j).name));
        for k=1:size(t.framewise_displacement,1)
            fd(end+1) = str2double(t.framewise_displacement(k,:));
        end
    end
    c.data(:,fd<0.3) = []; % remove high motion volumes

    % run pfm and filter
    Densities=flip([0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05]);
    NumberReps=[1 5 10 50 50 50 50 50 50];
    Structures={'CORTEX_LEFT','CEREBELLUM_LEFT','ACCUMBENS_LEFT','CAUDATE_LEFT','PALLIDUM_LEFT','PUTAMEN_LEFT','THALAMUS_LEFT','HIPPOCAMPUS_LEFT','AMYGDALA_LEFT','ACCUMBENS_LEFT',...
        'CORTEX_RIGHT','CEREBELLUM_RIGHT','ACCUMBENS_RIGHT','CAUDATE_RIGHT','PALLIDUM_RIGHT','PUTAMEN_RIGHT','THALAMUS_RIGHT','HIPPOCAMPUS_RIGHT','AMYGDALA_RIGHT','ACCUMBENS_RIGHT'};
    MinDistance = 10;
    BadVerts = [];
    NumberCores = {pc cores};
    pfm(c,distMat,tmpDir,Densities,NumberReps,MinDistance,BadVerts,Structures,NumberCores)
    spatial_filtering(fullfile(tmpDir,'Bipartite_PhysicalCommunities.dtseries.nii'),outDir,sprintf('sub-%s_ses-%s_Bipartite_PhysicalCommunities_desc-SpatialFiltering.dtseries.nii',sub,ses{i}),...
        surf,inputs.spatialThresh(1),inputs.spatialThresh(2));
end

end
