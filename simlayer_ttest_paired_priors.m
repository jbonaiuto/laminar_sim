function simlayer_ttest_paired_priors( subj_info, session_num, invfoi, SNR, varargin )

% Parse inputs
defaults = struct('surf_dir', 'd:\pred_coding\surf', 'mri_dir', 'd:\pred_coding\mri',...
    'dipole_moment', 10);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Copy already-inverted file
rawfile=fullfile('/data/pred_coding/analysis/',subj_info.subj_id, num2str(session_num), 'grey_coreg\EBB\p0.4\instr\f15_30', sprintf('r%s_%d.mat',subj_info.subj_id,session_num));
% Output directory
out_path=fullfile('/data/layer_sim/ttest_paired_prior_results',subj_info.subj_id,num2str(session_num),sprintf('f%d_%d_SNR%d_dipolemoment%d',invfoi(1),invfoi(2),SNR,params.dipole_moment));
if exist(out_path,'dir')~=7
    mkdir(out_path);
end
% New file to work with
newfile=fullfile(out_path, sprintf('%s_%d.mat',subj_info.subj_id,session_num));
greyregfile=fullfile(out_path, sprintf('%s_%d_greycoreg.mat',subj_info.subj_id,session_num));

spm('defaults', 'EEG');
spm_jobman('initcfg'); 

% Copy file to foi_dir
clear jobs
matlabbatch=[];
matlabbatch{1}.spm.meeg.other.copy.D = {rawfile};
matlabbatch{1}.spm.meeg.other.copy.outfile = newfile;
spm_jobman('run', matlabbatch);

% Copy file to foi_dir
clear jobs
matlabbatch=[];
matlabbatch{1}.spm.meeg.other.copy.D = {rawfile};
matlabbatch{1}.spm.meeg.other.copy.outfile = greyregfile;
spm_jobman('run', matlabbatch);

% Load meshes
orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');
orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');
pialwhite_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed-ds_pial.hires.deformed.surf.gii');
simmeshes={white_mesh,pial_mesh};
allmeshes={white_mesh,pial_mesh,pialwhite_mesh};

% Create smoothed meshes
patch_extent_mm=5; %5 approx mm
for meshind=1:length(allmeshes),
    [path,file,ext]=fileparts(allmeshes{meshind});
    smoothedfile=fullfile(path, sprintf('FWHM%d.00_%s.mat',patch_extent_mm,file));
    if exist(smoothedfile,'file')~=2
        tic
        [smoothkern]=spm_eeg_smoothmesh_mm(allmeshes{meshind},abs(patch_extent_mm));
        toc
    end
end
pial=gifti(pial_mesh);
white=gifti(white_mesh);

pial_white_map=map_pial_to_white(white_mesh, pial_mesh, ...
        'mapType', 'link', 'origPial', orig_pial_mesh, ...
        'origWhite', orig_white_mesh);
    
white_pial_map=map_white_to_pial(white_mesh, pial_mesh, ...
        'mapType', 'link', 'origPial', orig_pial_mesh, ...
        'origWhite', orig_white_mesh);

%% Setup simulation - number of sources, list of vertices to simulate on
nverts=size(white.vertices,1);
rng(0);
simvertind=randperm(nverts); %% random list of vertex indices to simulate sources on
Nsim=60; %% number of simulated sources on each surface

%% for MSP  or GS or ARD
% Number of patches as priors
% so use all vertices that will be simulated on, on white/pial surface
% (plus a few more) as MSP priors
Ip=[simvertind(1:Nsim) nverts+white_pial_map(simvertind(1:Nsim))' pial_white_map(simvertind(1:Nsim))' nverts+simvertind(1:Nsim)];
% Save priors
patchfilename=fullfile(out_path, 'temppatch.mat');
save(patchfilename,'Ip');

methodnames={'MSP'};
Nmeth=length(methodnames);

% Inversion parameters
invwoi=[-500 500];
simwoi=[100 500];
baselinewoi=[-500 -100];
% Number of cross validation folds
Nfolds=1;
% Percentage of test channels in cross validation
ideal_pctest=0;
% Use all available spatial modes
ideal_Nmodes=[];

% Coregister simulated dataset to combined pial/white mesh
matlabbatch=[];
matlabbatch{1}.spm.meeg.source.headmodel.D = {greyregfile};
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.mri_dir,[subj_info.subj_id subj_info.birth_date], [subj_info.headcast_t1 ',1'])};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {pialwhite_mesh};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = subj_info.nas;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = subj_info.lpa;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = subj_info.rpa;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
spm_jobman('run',matlabbatch);
        
% Setup spatial modes for cross validation
spatialmodesname=fullfile(out_path, 'testmodes.mat');
[spatialmodesname,Nmodes,pctest]=spm_eeg_inv_prep_modes_xval(greyregfile, ideal_Nmodes, spatialmodesname, Nfolds, ideal_pctest);

greycoreg=load(greyregfile);

for simmeshind=1:length(simmeshes)
    simmesh=simmeshes{simmeshind};

    % coregister to correct mesh
    spm_jobman('initcfg'); 
    matlabbatch=[];
    matlabbatch{1}.spm.meeg.source.headmodel.D = {newfile};
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.mri_dir,[subj_info.subj_id subj_info.birth_date], [subj_info.headcast_t1 ',1'])};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {simmesh};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = subj_info.nas;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = subj_info.lpa;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = subj_info.rpa;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    spm_jobman('run',matlabbatch);
        
    Dmesh=spm_eeg_load(newfile);
    
    %% now simulate sources on this mesh
    for s=1:Nsim,
        %% get location to simulate dipole on this mesh
        simpos=Dmesh.inv{1}.mesh.tess_mni.vert(simvertind(s),:); 

        prefix=sprintf('sim_mesh%d_source%d',simmeshind, s);

        % Simulate source 
        matlabbatch=[];
        matlabbatch{1}.spm.meeg.source.simulate.D = {newfile};
        matlabbatch{1}.spm.meeg.source.simulate.val = 1;
        matlabbatch{1}.spm.meeg.source.simulate.prefix = prefix;
        matlabbatch{1}.spm.meeg.source.simulate.whatconditions.all = 1;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.woi = simwoi;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.isSin.foi = mean(invfoi);
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.dipmom = [params.dipole_moment patch_extent_mm];
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.locs = simpos;
        if abs(params.dipole_moment)>0
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.setSNR = SNR;               
        else
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.whitenoise = 100;
        end
        [a,b]=spm_jobman('run', matlabbatch);
        
        % Copy forward model from pial/white coregistered file
        simfilename=fullfile(out_path,sprintf('%s%s_%d.mat',prefix,subj_info.subj_id,session_num));
        sim=load(simfilename);
        sim.D.other=greycoreg.D.other;
        D=sim.D;
        copyfile(fullfile(out_path, sprintf('SPMgainmatrix_%s_%d_greycoreg_1.mat', subj_info.subj_id, session_num)), fullfile(out_path, sprintf('SPMgainmatrix_%s%s_%d_1.mat', prefix, subj_info.subj_id, session_num)));
        D.other.inv{1}.gainmat=sprintf('SPMgainmatrix_%s%s_%d_1.mat', prefix, subj_info.subj_id, session_num);
        save(simfilename,'D');
                    
        % Resconstruct using each method
        for methind=1:Nmeth,       
            method=methodnames{methind};
            
            % Do inversion of simulated data with this surface
            matlabbatch=[];
            matlabbatch{1}.spm.meeg.source.invertiter.D = {simfilename};
            matlabbatch{1}.spm.meeg.source.invertiter.val = 1;
            matlabbatch{1}.spm.meeg.source.invertiter.whatconditions.all = 1;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invtype = method;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.woi = invwoi;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.foi = invfoi;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.hanning = 0;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedfile = {patchfilename}; % '<UNDEFINED>';
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedrows = 1; %'<UNDEFINED>';
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm =[-patch_extent_mm]; %% NB A fiddle here- need to properly quantify
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = Nmodes;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.umodes = {spatialmodesname};
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = [];
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
            matlabbatch{1}.spm.meeg.source.invertiter.modality = {'All'};
            matlabbatch{1}.spm.meeg.source.invertiter.crossval = [pctest Nfolds];   
            spm_jobman('run',matlabbatch);

            matlabbatch=[];
            % Extract results and smooth
            matlabbatch{1}.spm.meeg.source.results.D = {simfilename};
            matlabbatch{1}.spm.meeg.source.results.val = 1;
            matlabbatch{1}.spm.meeg.source.results.woi = [baselinewoi; simwoi];
            matlabbatch{1}.spm.meeg.source.results.foi = invfoi;
            matlabbatch{1}.spm.meeg.source.results.ctype = 'trials';
            matlabbatch{1}.spm.meeg.source.results.space = 0;
            matlabbatch{1}.spm.meeg.source.results.format = 'mesh';
            matlabbatch{1}.spm.meeg.source.results.smoothing = 8;
            spm_jobman('run',matlabbatch);

            % Split sources in WOI
            woi_dir=fullfile(out_path, ['t' num2str(simwoi(1)) '_' num2str(simwoi(2))]);
            mkdir(woi_dir);
            delete(fullfile(woi_dir,'*'));
            movefile(fullfile(out_path, sprintf('%s%s_%d_1_t%d_%d_f%d_%d_*', prefix, subj_info.subj_id, session_num, simwoi(1), simwoi(2), invfoi(1), invfoi(2))), woi_dir);
            split_inversion_results(woi_dir);

            % Split sources in baseline
            baseline_dir=fullfile(out_path, ['t' num2str(baselinewoi(1)) '_' num2str(baselinewoi(2))]);
            mkdir(baseline_dir);
            delete(fullfile(baseline_dir,'*'));
            movefile(fullfile(out_path, sprintf('%s%s_%d_1_t%d_%d_f%d_%d_*', prefix, subj_info.subj_id, session_num, baselinewoi(1), baselinewoi(2), invfoi(1), invfoi(2))), baseline_dir);
            split_inversion_results(baseline_dir);

            % Load all pial data from wois
            pial_prefix=sprintf('pial_%s%s_%d_1_', prefix, subj_info.subj_id, session_num);
            pial_woi_trials=load_woi_trials(woi_dir, pial_prefix, simwoi, invfoi, size(pial.vertices,1));
            pial_baseline_trials=load_woi_trials(baseline_dir, pial_prefix, baselinewoi, invfoi, size(pial.vertices,1));

            % Load all white matter data from wois
            white_prefix=sprintf('white_%s%s_%d_1_', prefix, subj_info.subj_id, session_num);
            white_woi_trials=load_woi_trials(woi_dir, white_prefix, simwoi, invfoi, size(white.vertices,1));
            white_baseline_trials=load_woi_trials(baseline_dir, white_prefix, baselinewoi, invfoi, size(white.vertices,1));

            % Save pial diff
            pial_diff=pial_woi_trials-pial_baseline_trials;
            white_diff=white_woi_trials-white_baseline_trials;

            write_metric_gifti(fullfile(out_path, sprintf('pial.%s.%s.gii',method,prefix)), pial_diff);
            write_metric_gifti(fullfile(out_path, sprintf('white.%s.%s.gii',method,prefix)), white_diff);

        end
        close all
        delete(fullfile(out_path, sprintf('sim_mesh%d_source%d*.*',simmeshind, s)));
    end
end
