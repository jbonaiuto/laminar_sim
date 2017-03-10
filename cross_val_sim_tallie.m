function cross_val_sim_tallie(subj_info, session_num, invfoi, SNR, varargin)

% Parse inputs
defaults = struct('surf_dir', 'd:\pred_coding\surf', 'mri_dir', 'd:\pred_coding\mri',...
    'out_file', '', 'dipole_moment', 10);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Copy already-inverted file
rawfile=fullfile('D:\pred_coding\analysis\',subj_info.subj_id, num2str(session_num), 'grey_coreg\EBB\p0.4\instr\f15_30', sprintf('r%s_%d.mat',subj_info.subj_id,session_num));
% Output directory
out_path=fullfile('D:\layer_sim\results',subj_info.subj_id,num2str(session_num));
if exist(out_path,'dir')~=7
    mkdir(out_path);
end
% New file to work with
newfile=fullfile(out_path, sprintf('%s_%d.mat',subj_info.subj_id,session_num));

spm('defaults', 'EEG');
spm_jobman('initcfg'); 

% Copy file to foi_dir
clear jobs
matlabbatch=[];
matlabbatch{1}.spm.meeg.other.copy.D = {rawfile};
matlabbatch{1}.spm.meeg.other.copy.outfile = newfile;
spm_jobman('run', matlabbatch);

% meanpos=mean(M.vertices);
% newvert=M.vertices-repmat(meanpos,length(M.vertices),1);
% newvert=newvert./2+repmat(meanpos,length(M.vertices),1);
% M.vertices=newvert;
% save(M,'D:\matlab\batch\halfmesh.gii');


% White and pial meshes for this subject
allmeshes=strvcat(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii'),...
    fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii'));


Nmesh=size(allmeshes,1);

% Create smoothed meshes
patch_extent_mm=-5; %5 approx mm
for meshind=1:Nmesh,
    [path,file,ext]=fileparts(deblank(allmeshes(meshind,:)));
    smoothedfile=fullfile(path, sprintf('FWHM5.00_%s.mat',file));
    if exist(smoothedfile,'file')~=2
        tic
        [smoothkern]=spm_eeg_smoothmesh_mm(deblank(allmeshes(meshind,:)),abs(patch_extent_mm));
        toc
    end
end

%% Setup simulation - number of sources, list of vertices to simulate on
mesh_one=gifti(allmeshes(1,:));
nverts=size(mesh_one.vertices,1);
rng(0);
simvertind=randperm(nverts); %% random list of vertex indices to simulate sources on
%Nsim=8; %% number of simulated sources
Nsim=60; %% number of simulated sources

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
Nfolds=10;
% Percentage of test channels in cross validation
ideal_pctest=10; %% may not use this number as we need integer number of channels
% Use all available spatial modes
ideal_Nmodes=[];


% All F values and cross validation errors
% meshes simulated on x number of simulations x meshes reconstructed onto x
% num methods x num cross validation folds
allcrossErr=zeros(Nmesh,Nsim,Nmesh,Nmeth,Nfolds);


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
        
        % Simulate source 
        matlabbatch=[];
        matlabbatch{1}.spm.meeg.source.simulate.D = {filename};
        matlabbatch{1}.spm.meeg.source.simulate.val = 1;
        matlabbatch{1}.spm.meeg.source.simulate.prefix = sprintf('sim_mesh%d_source%d',simmeshind,s);
        matlabbatch{1}.spm.meeg.source.simulate.whatconditions.all = 1;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.woi = invwoi;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.isSin.foi = mean(invfoi);
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.dipmom = [params.dipole_moment patch_extent_mm];
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.locs = simpos;
        if abs(params.dipole_moment)>0
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.setSNR = SNR;               
        else
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.whitenoise = 100;
        end
        [a,b]=spm_jobman('run', matlabbatch);
        
        % Load simulated dataset
        simfilename=a{1}.D{1};        
        Dsim=spm_eeg_load(simfilename);        
        
        %% now reconstruct onto all the meshes and look at cross val and F vals
        for meshind=1:Nmesh,
            
            % Coregister simulated dataset to reconstruction mesh
            matlabbatch=[];
            matlabbatch{1}.spm.meeg.source.headmodel.D = {simfilename};
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
            
            % Setup spatial modes for cross validation
            spatialmodesname=[Dsim.path filesep 'testmodes.mat'];
            [spatialmodesname,Nmodes,pctest]=spm_eeg_inv_prep_modes_xval(simfilename, ideal_Nmodes, spatialmodesname, Nfolds, ideal_pctest);
            
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
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm =[patch_extent_mm]; %% NB A fiddle here- need to properly quantify
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
                [a1,b1]=spm_jobman('run', matlabbatch);
                
                % Load inversion - get cross validation error end F
                Drecon=spm_eeg_load(simfilename);                
                allcrossErr(simmeshind,s,meshind,methind,:)=Drecon.inv{1}.inverse.crosserr;
                
                
            end; % for methind                        
        end; %% for reconstruction mesh (meshind)
        close all;
    end; % for s (sources)
end; % for simmeshind (simulatiom mesh)
if length(params.out_file)==0
    params.out_file=sprintf('allcrossErr_f%d_%d_SNR%d.mat',invfoi(1),invfoi(2),SNR);
end
save(fullfile(out_path,params.out_file),'allcrossErr');

for methind=1:Nmeth,                
    figure(methind);clf;
%   figure(2);clf;

    % For each simulated mesh
    for simmeshind=1:Nmesh,
        [path,file,ext]=fileparts(deblank(allmeshes(simmeshind,:)));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        simmeshname=y{2};
        % other mesh index (assuming there are just 2 meshes)
        otherind=setxor(simmeshind,1:Nmesh);
        [path,file,ext]=fileparts(deblank(allmeshes(otherind,:)));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        othermeshname=y{2};

        % F reconstructed on true - reconstructed on other
        % num simulations x number of folds
        %truotherF=squeeze(allcrossF(simmeshind,:,simmeshind,methind,:)-allcrossF(simmeshind,:,otherind,methind,:));
        truotherF=squeeze(mean(allcrossErr(simmeshind,:,otherind,methind,:),5)-mean(allcrossErr(simmeshind,:,simmeshind,methind,:),5));
        %figure(1);
        subplot(Nmesh,1,simmeshind);
        bar(truotherF)
        xlabel('Simulation')
        ylabel('Crossval Err Difference');
        title(sprintf('Crossval Error, %s, %s-%s',methodnames{methind},simmeshname,othermeshname));        

        % Cross val reconstructed on other -reconstructed on true ( should be positive)
%         truotherXval=squeeze(allcrosserr(simmeshind,:,otherind,methind,:)-allcrosserr(simmeshind,:,simmeshind,methind,:));
%         figure(2);
%         subplot(Nmesh,1,simmeshind);    
%         bar(truotherXval);
%         xlabel('Simulation');
%         ylabel('X val error diff');
%         title(sprintf('Cross val, %s-%s',othermeshname,simmeshname));
%         % Average over simulations and folds
%         meanXval(simmeshind)=mean(mean(truotherXval));    
%         % Number of correct cross validations (err on other > err on true)
%         correctXval(simmeshind)=length(find(truotherXval>0))/length(truotherXval(:));
    end

end
% for methind=1:Nmeth,                
%     truotherF=squeeze(allcrossF(simmeshind,:,simmeshind,methind)-allcrossF(simmeshind,:,otherind,methind));
%     % Average over folds
%     %meanF((simmeshind-1)*Nsim+1:simmeshind*Nsim)=mean(truotherF,2)';        
%     
%     %% do the stats on the group
%     lme=[zeros(Nsim*Nmesh,1) truotherF'] %% difference between true and other mesh
%     spm_BMS(lme,[],1)
% end


%error('stop');

% meshes simulated on x number of simulations x meshes reconstructed onto x
% num methods x num cross validation folds

% figure;
% plot(squeeze(allbatF(:,1,:)),squeeze(allbaterr(:,1,:)),'x');
% hold on;
% plot(squeeze(allbatF(:,2,:)),squeeze(allbaterr(:,2,:)),'o');
% 
% 
% 
% 
% crossval=zeros(Nmodels,Nmodels);
% 
% Fdiff=zeros(Nmodels,Nmodels);
% Flmisdiff=[];
% for j=1:Nmodels,
%     for k=1:Nmodels,
%         Fdiff(j,k)=mean(allbatF(j,1,:)-allbatF(k,1,:));
%         
%         
%         errdiff=allbaterr(j,1,:)-allbaterr(k,1,:);
%         
%         crossval(j,k)=length(find(errdiff<0))./Ntest; %% where j is better than k
%         
%         
%     end; %for k
% end; % for j
% figure;
% subplot(3,1,1);
% imagesc(Fdiff);colorbar;
% subplot(3,1,2)
% imagesc(crossval);colorbar;
% 
% 
% dum=triu(ones(size(crossval)));
% flatind=find(dum);
% figure;
% 
% 
% probF =1./(1 + exp(-Fdiff(flatind)));  %% 300ms of data
% 
% 
% figure;
% plot(Fdiff(flatind),crossval(flatind),'g.');% 'r.',Fdiff(flatind),probF,'g.');
% legend('cross val','F pval');
% 
% figure;
% plot(Fdiff(flatind),crossval2(flatind),'r.'); %% (l % Fdiff(flatind),probF,'g.');
% legend('cross val 2','F pval');
% 
% 
% 
