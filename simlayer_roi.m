function simlayer_roi( subj_info, session_num, invfoi, SNR, varargin )
% SIMLAYER_ROI  Run simulations with ROI analysis
%
% Use as
%   simlayer_roi(subjects(1), 1, [10 30], -20)
% where the first argument is the subject info structure (from create_subjects),
% the second is the session number, the third is the frequency range, and
% the fourth is the SNR (db).
% 
%   simlayer_roi(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * surf_dir - directory containing subject surfaces
%    * mri_dir - directory containing subject MRIs
%    * out_path - output file path (automatically generated if not
%    specified)
%    * dipole_moment - 10 (default) or interger - moment of simulated
%    dipole
%    * sim_patch_size - 5 (default) or interger - simulated patch size
%    * reconstruct_patch_size - 5 (default) or interger - reconstruction patch size

% Parse inputs
defaults = struct('surf_dir', 'd:\pred_coding\surf', 'mri_dir', 'd:\pred_coding\mri',...
    'out_path', '', 'dipole_moment', 10, 'sim_patch_size', 5,...
    'reconstruct_patch_size', 5);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Copy already-inverted file
rawfile=fullfile('d:/pred_coding/analysis',subj_info.subj_id,...
    num2str(session_num), 'grey_coreg\EBB\p0.4\instr\f15_30',...
    sprintf('br%s_%d.mat',subj_info.subj_id,session_num));
% Output directory
if length(params.out_path)==0
    params.out_path=fullfile('d:/layer_sim/ttest_results',subj_info.subj_id,...
        num2str(session_num),sprintf('f%d_%d_SNR%d_dipolemoment%d',invfoi(1),...
        invfoi(2),SNR,params.dipole_moment));
end
if exist(params.out_path,'dir')~=7
    mkdir(params.out_path);
end
% New file to work with
newfile=fullfile(params.out_path, sprintf('%s_%d.mat',subj_info.subj_id,session_num));
greyregfile=fullfile(params.out_path, sprintf('%s_%d_greycoreg.mat',subj_info.subj_id,session_num));

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
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');
pialwhite_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed-ds_pial.hires.deformed.surf.gii');
simmeshes={white_mesh,pial_mesh};
allmeshes={white_mesh,pial_mesh,pialwhite_mesh};

% Create smoothed meshes
for meshind=1:length(allmeshes),
    [smoothkern]=spm_eeg_smoothmesh_mm(allmeshes{meshind},params.sim_patch_size);
end
pial=gifti(pial_mesh);
white=gifti(white_mesh);

%% Setup simulation - number of sources, list of vertices to simulate on
nverts=size(white.vertices,1);
rng(0);
simvertind=randperm(nverts); %% random list of vertex indices to simulate sources on
Nsim=60; %% number of simulated sources on each surface

%% for MSP  or GS or ARD
% Number of patches as priors
Npatch=Nsim*1.25;
% so use all vertices that will be simulated on, on white/pial surface
% (plus a few more) as MSP priors
Ip=[simvertind(1:Npatch) nverts+simvertind(1:Npatch)];
% Save priors
patchfilename=fullfile(params.out_path, 'temppatch.mat');
save(patchfilename,'Ip');

methodnames={'EBB','IID','COH','MSP'};
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
spatialmodesname=fullfile(params.out_path, 'testmodes.mat');
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
    sims=[1:Nsim];
    for s=sims,
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
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.dipmom = [params.dipole_moment params.sim_patch_size];
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.locs = simpos;
        if abs(params.dipole_moment)>0
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.setSNR = SNR;               
        else
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.whitenoise = 100;
        end
        [a,b]=spm_jobman('run', matlabbatch);
        
        % Copy forward model from pial/white coregistered file
        simfilename=fullfile(params.out_path,sprintf('%s%s_%d.mat',prefix,subj_info.subj_id,session_num));
        sim=load(simfilename);
        sim.D.other=greycoreg.D.other;
        D=sim.D;
        copyfile(fullfile(params.out_path, sprintf('SPMgainmatrix_%s_%d_greycoreg_1.mat', subj_info.subj_id, session_num)), fullfile(params.out_path, sprintf('SPMgainmatrix_%s%s_%d_1.mat', prefix, subj_info.subj_id, session_num)));
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
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm =[-params.reconstruct_patch_size]; %% NB A fiddle here- need to properly quantify
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

            D=spm_eeg_load(simfilename);
            goodchans=D.indchantype('MEGGRAD','good');
            M=D.inv{1}.inverse.M;
            U=D.inv{1}.inverse.U{1};
            T=D.inv{1}.inverse.T;
            It   = D.inv{1}.inverse.It;
            times=D.inv{1}.inverse.pst;
            Dgood=squeeze(D(goodchans,It,:));
            ntrials=size(Dgood,3);
            m = export(gifti(D.inv{1}.mesh.tess_ctx),'patch');
            GL      = spm_mesh_smooth(m);

            wois=[baselinewoi; simwoi];
            for w=1:size(wois,1)
                woi=wois(w,:);

                fwhm = max(diff(woi),8);
                t    = exp(-4*log(2)*(times(:) - mean(woi)).^2/(fwhm^2));
                t    = t/sum(t);

                % get frequency space and put PST subspace into contrast (W -> T*T'*W)
                %--------------------------------------------------------------------------
                wt = 2*pi*times(:)/1000;
                W  = [];
                for f = invfoi(1):invfoi(end)
                    W = [W sin(f*wt) cos(f*wt)];
                end
                W  = diag(t)*W;
                W  = spm_svd(W,1);  
                TW     = T'*W;
                TTW    = T*TW;

                woi_vals=zeros(nverts*2,ntrials);
                MU=M*U;
                for i=1:ntrials
                    MUd1=MU*squeeze(Dgood(:,:,i));
                    Y     = sum((MUd1*TTW).^2,2);
                    woi_vals(:,i)=spm_mesh_smooth(GL,Y,8);       
                end

                woi_dir=fullfile(params.out_path, ['t' num2str(woi(1)) '_' num2str(woi(2))]);
                if exist(woi_dir,'dir')~=7
                    mkdir(woi_dir);
                end
                delete(fullfile(woi_dir,'*'));
                out_filename=fullfile(woi_dir, sprintf('%s%s_%d_1_t%d_%d_f%d_%d', prefix, subj_info.subj_id, session_num, woi(1), woi(2), invfoi(1), invfoi(2)));
                write_metric_gifti(out_filename,woi_vals);    

                % Split pial and grey sources
                split_inversion_results(woi_dir);

            end
            
            baseline_dir=fullfile(params.out_path, ['t' num2str(baselinewoi(1)) '_' num2str(baselinewoi(2))]);
            woi_dir=fullfile(params.out_path, ['t' num2str(simwoi(1)) '_' num2str(simwoi(2))]);
            
            % Load all pial data from wois
            pial_woi_trials=gifti(fullfile(woi_dir,sprintf('pial_%s%s_%d_1_t%d_%d_f%d_%d.gii', prefix, subj_info.subj_id, session_num, simwoi(1), simwoi(2), invfoi(1), invfoi(2))));
            pial_baseline_trials=gifti(fullfile(baseline_dir,sprintf('pial_%s%s_%d_1_t%d_%d_f%d_%d.gii', prefix, subj_info.subj_id, session_num, baselinewoi(1), baselinewoi(2), invfoi(1), invfoi(2))));

            % Load all white matter data from wois
            white_woi_trials=gifti(fullfile(woi_dir,sprintf('white_%s%s_%d_1_t%d_%d_f%d_%d.gii', prefix, subj_info.subj_id, session_num, simwoi(1), simwoi(2), invfoi(1), invfoi(2))));
            white_baseline_trials=gifti(fullfile(baseline_dir,sprintf('white_%s%s_%d_1_t%d_%d_f%d_%d.gii', prefix, subj_info.subj_id, session_num, baselinewoi(1), baselinewoi(2), invfoi(1), invfoi(2))));

            % Save pial diff
            pial_diff=pial_woi_trials.cdata(:,:)-pial_baseline_trials.cdata(:,:);
            white_diff=white_woi_trials.cdata(:,:)-white_baseline_trials.cdata(:,:);

            write_metric_gifti(fullfile(params.out_path, sprintf('pial.%s.%s.gii',method,prefix)), pial_diff);
            write_metric_gifti(fullfile(params.out_path, sprintf('white.%s.%s.gii',method,prefix)), white_diff);

        end
        close all
        delete(fullfile(params.out_path, sprintf('sim_mesh%d_source%d*.*',simmeshind, s)));
    end
end
