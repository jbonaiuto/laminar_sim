function simlayer_free_energy(subj_info, session_num, invfoi, SNR, varargin)
% SIMLAYER_FREE_ENERGY  Run simulations with whole brain - free energy
%   analysis
%
% Use as
%   simlayer_free_energy(subjects(1), 1, [10 30], -20)
% where the first argument is the subject info structure (from create_subjects),
% the second is the session number, the third is the frequency range, and
% the fourth is the SNR (db).
% 
%   simlayer_free_energy(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * surf_dir - directory containing subject surfaces
%    * mri_dir - directory containing subject MRIs
%    * out_file - output file name (automatically generated if not
%    specified)
%    * dipole_moment - 10 (default) or interger - moment of simulated
%    dipole
%    * sim_patch_size - 5 (default) or interger - simulated patch size
%    * reconstruct_patch_size - 5 (default) or interger - reconstruction patch size
%    * nsims - 60 (default) or integer - number of simulations per surface

% Parse inputs
defaults = struct('surf_dir', 'd:\pred_coding\surf', 'mri_dir', 'd:\pred_coding\mri',...
    'out_file', '', 'dipole_moment', 10, 'sim_patch_size', 5,...
    'reconstruct_patch_size', 5, 'nsims', 60);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Copy already-inverted file
rawfile=fullfile('d:/pred_coding/analysis/',subj_info.subj_id,...
    num2str(session_num), 'grey_coreg\EBB\p0.4\instr\f15_30',...
    sprintf('br%s_%d.mat',subj_info.subj_id,session_num));
% Output directory
out_path=fullfile('c:/layer_sim/results',subj_info.subj_id,num2str(session_num));
if exist(out_path,'dir')~=7
    mkdir(out_path);
end
% New file to work with
newfile=fullfile(out_path, sprintf('%s_%d.mat',subj_info.subj_id,session_num));

if length(params.out_file)==0
    params.out_file=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',...
        invfoi(1),invfoi(2),SNR,params.dipole_moment);
end

spm('defaults', 'EEG');
spm_jobman('initcfg'); 

% Copy file to foi_dir
clear jobs
matlabbatch=[];
matlabbatch{1}.spm.meeg.other.copy.D = {rawfile};
matlabbatch{1}.spm.meeg.other.copy.outfile = newfile;
spm_jobman('run', matlabbatch);

% White and pial meshes for this subject
allmeshes=strvcat(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii'),...
    fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii'));
Nmesh=size(allmeshes,1);

% Create smoothed meshes
for meshind=1:Nmesh,
    [smoothkern]=spm_eeg_smoothmesh_mm(deblank(allmeshes(meshind,:)),params.sim_patch_size);
end

%% Setup simulation - number of sources, list of vertices to simulate on
mesh_one=gifti(allmeshes(1,:));
nverts=size(mesh_one.vertices,1);
rng(0);
simvertind=randperm(nverts); %% random list of vertex indices to simulate sources on
Nsim=params.nsims; %% number of simulated sources

%% for MSP  or GS or ARD
% Number of patches as priors
Npatch=round(Nsim*1.5);
% so use all vertices that will be simulated on (plus a few more) as MSP priors
Ip=simvertind(1:Npatch);

% Save priors
patchfilename=fullfile(out_path, 'temppatch.mat');
save(patchfilename,'Ip');

% Inversion method to use
methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
Nmeth=length(methodnames);

% Inversion parameters
invwoi=[100 500];
% Number of cross validation folds
Nfolds=1;
% Percentage of test channels in cross validation
ideal_pctest=0;
% Use all available spatial modes
ideal_Nmodes=[];


% All F values
allcrossF=zeros(Nmesh,Nsim,Nmesh,Nmeth);
allcrossVE=zeros(Nmesh,Nsim,Nmesh,Nmeth);

regfiles={};
spatialmodesnames={};

for meshind=1:Nmesh,
    regfile=fullfile(out_path, sprintf('%s_%d_%dcoreg.mat',subj_info.subj_id,session_num,meshind));
    regfiles{meshind}=regfile;
    if exist(regfile,'file')~=2
        clear jobs
        matlabbatch=[];
        matlabbatch{1}.spm.meeg.other.copy.D = {rawfile};
        matlabbatch{1}.spm.meeg.other.copy.outfile = regfile;
        spm_jobman('run', matlabbatch);
        
        % Coregister simulated dataset to reconstruction mesh
        matlabbatch=[];
        matlabbatch{1}.spm.meeg.source.headmodel.D = {regfile};
        matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
        matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.mri_dir,[subj_info.subj_id subj_info.birth_date], [subj_info.headcast_t1 ',1'])};
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {deblank(allmeshes(meshind,:))};
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
        spm_jobman('run', matlabbatch);                   
    end
    % Setup spatial modes for cross validation
    spatialmodesname=fullfile(out_path, sprintf('%d_testmodes.mat',meshind));    
    [spatialmodesname,Nmodes,pctest]=spm_eeg_inv_prep_modes_xval(regfile, ideal_Nmodes, spatialmodesname, Nfolds, ideal_pctest);
    spatialmodesnames{meshind}=spatialmodesname;
end


% Simulate sources on each mesh
for simmeshind=1:Nmesh, %% choose mesh to simulate on
    
    simmesh=deblank(allmeshes(simmeshind,:));
    
    %% coregister to correct mesh
    filename=deblank(newfile);
    matlabbatch=[];
    matlabbatch{1}.spm.meeg.source.headmodel.D = {filename};
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
    spm_jobman('run', matlabbatch);
    
    Dmesh=spm_eeg_load(filename);    
    
    %% now simulate sources on this mesh
    for s=1:Nsim,
        %% get location to simulate dipole on this mesh
        simpos=Dmesh.inv{1}.mesh.tess_mni.vert(simvertind(s),:); 
        prefix=sprintf('sim_mesh%d_source%d',simmeshind,s);

        % Simulate source 
        matlabbatch=[];
        matlabbatch{1}.spm.meeg.source.simulate.D = {filename};
        matlabbatch{1}.spm.meeg.source.simulate.val = 1;
        matlabbatch{1}.spm.meeg.source.simulate.prefix = prefix;
        matlabbatch{1}.spm.meeg.source.simulate.whatconditions.all = 1;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.woi = invwoi;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.isSin.foi = mean(invfoi);
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.dipmom = [params.dipole_moment params.sim_patch_size];
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.locs = simpos;
        if abs(params.dipole_moment)>0
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.setSNR = SNR;               
        else
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.whitenoise = 100;
        end
        [a,b]=spm_jobman('run', matlabbatch);
        
        %% now reconstruct onto all the meshes and look at cross val and F vals
        for meshind=1:Nmesh,
        
            % Copy forward model from pial or white coregistered file
            simfilename=fullfile(out_path,sprintf('%s%s_%d.mat',prefix,subj_info.subj_id,session_num));
            sim=load(simfilename);
            reconcoreg=load(regfiles{meshind});
            sim.D.other=reconcoreg.D.other;
            D=sim.D;
            copyfile(fullfile(out_path, sprintf('SPMgainmatrix_%s_%d_%dcoreg_1.mat', subj_info.subj_id, session_num, meshind)), fullfile(out_path, sprintf('SPMgainmatrix_%s%s_%d_1.mat', prefix, subj_info.subj_id, session_num)));
            D.other.inv{1}.gainmat=sprintf('SPMgainmatrix_%s%s_%d_1.mat', prefix, subj_info.subj_id, session_num);
            save(simfilename,'D');
    
            % Resconstruct using each method
            for methind=1:Nmeth,                
                
                % Do inversion of simulated data with this surface
                matlabbatch=[];
                matlabbatch{1}.spm.meeg.source.invertiter.D = {simfilename};
                matlabbatch{1}.spm.meeg.source.invertiter.val = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.whatconditions.all = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invtype = methodnames{methind}; %;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.woi = invwoi;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.foi = invfoi;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.hanning = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedfile = {patchfilename}; % '<UNDEFINED>';
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedrows = 1; %'<UNDEFINED>';
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm =[-params.reconstruct_patch_size]; %% NB A fiddle here- need to properly quantify
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = Nmodes;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.umodes = {spatialmodesnames{simmeshind}};
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = [];
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
                matlabbatch{1}.spm.meeg.source.invertiter.modality = {'All'};
                matlabbatch{1}.spm.meeg.source.invertiter.crossval = [pctest Nfolds];                                
                [a1,b1]=spm_jobman('run', matlabbatch);
                
                % Load inversion - get cross validation error end F
                Drecon=spm_eeg_load(simfilename);                
                allcrossF(simmeshind,s,meshind,methind)=Drecon.inv{1}.inverse.crossF;
                allcrossVE(simmeshind,s,meshind,methind)=Drecon.inv{1}.inverse.VE;
                
                
            end
        end
        close all;
    end
end
save(fullfile(out_path,params.out_file),'allcrossF','allcrossVE');


