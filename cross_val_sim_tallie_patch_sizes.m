% List of open inputs


inputs = cell(0, 0);
% TA % addpath D:\spm12

spm_jobman('initcfg')
spm('defaults', 'EEG');


% TA % rawfile= 'C:\Users\gbarnes\Documents\crossval\bemgfbinfspmeeg_gb070167_LUCIA_20120611_01.mat'; %% this can be any MEG dataset- it just gives head position, trial length , chan positions etc
rawfile = 'E:\MEG\s01\meffdspmeeg_SB260275_SheenaWaters_20150422_01_prepro_wp_f1-248_p0p4IID.mat';

% TA % M=gifti('D:\MEGgroup\freesurf_grbnative\surf\spm_lhrh_white.gii');
M = gifti('E:\MEG\s01\Fs_results\surf\ds_white.hires.deformed.surf.gii');
% meanpos=mean(M.vertices);
% newvert=M.vertices-repmat(meanpos,length(M.vertices),1);
% newvert=newvert./2+repmat(meanpos,length(M.vertices),1);
% M.vertices=newvert;
% save(M,'D:\matlab\batch\halfmesh.gii');


%allmeshes=strvcat('D:\MEGgroup\freesurf_grbnative\surf\spm_lhrh_white.gii',...
%    'D:\matlab\batch\halfmesh.gii'); %% white and pial meshes for this subject

% TA % allmeshes=strvcat('D:\MEGgroup\freesurf_grbnative\surf\spm_lhrh_white.gii',...
% TA %     'D:\matlab\batch\halfmesh.gii');
    %'D:\MEGgroup\freesurf_grbnative\surf\spm_lhrh_pial.gii');
allmeshes = strvcat('E:\MEG\s01\Fs_results\surf\ds_white.hires.deformed.surf.gii',...
    'E:\MEG\s01\Fs_results\surf\ds_pial.hires.deformed.surf.gii');


Nmesh=size(allmeshes,1);


D=spm_eeg_load(deblank(rawfile));
vert=D.inv{1}.mesh.tess_mni.vert;
rng(0);
simvertind=randperm(size(vert,1)/2); % TA "/2" %% random list of vertex indices to simulate sources on
Nsim=8; %% number of simulated sources

Npatch=30; %% for MSP  or GS or ARD
Ip=simvertind(1:Npatch); %% so use all vertices that will be simulated on (plus a few more) as MSP priors
patchfilename=[D.path filesep 'temppatch.mat'];
save(patchfilename,'Ip');


methodnames=strvcat('EBB'); %% just 1 method for now
Nmeth=size(methodnames,1);

invwoi=[100 500];
invfoi=[0 80];
SNR=5; %% dB
Nblocks=4;
ideal_pctest=10; %% may not use this number as we need integer number of channels
ideal_Nmodes=[];

patch_extent=0.75; %5 approx mm

allcrossF=zeros(Nmesh,Nsim,Nmesh,Nmeth,Nblocks);
allcrosserr=allcrossF;


for simmeshind=1:Nmesh, %% choose mesh to simulate on
    
    
    simmesh=deblank(allmeshes(simmeshind,:));
    
    %% coregister to correct mesh
    filename=deblank(rawfile);
    
    
    matlabbatch=[];
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.D = {filename};
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.comment = '';
    % TA % matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.mri = {'D:\MEGgroup\luzia\MQ0484.4\processed\msMQ0484-0004-00001-000208-01.img,1'};
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.mri = {'E:\MEG\s01\structural\s2014-07-11_11-51-115426-00001-00192-1_STRUCT.nii'};
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.cortex = {simmesh};
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.iskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.oskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.scalp = {''};
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshres = 2;
    % TA ? % matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregdefault = {'D:\MEGgroup\luzia\MQ0484.4\processed\helmet_coreg3.mat'};
    % TA:
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).fidname = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).specification.type = [-2.0190 110.5676 47.4948];
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).specification.type = [-83.6890 36.2986 -18.7222];
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).specification.type = [79.3990 40.6376 -17.9802];
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.useheadshape = 0;
    %
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.eeg = 'EEG BEM';
    matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.meg = 'Single Shell';
    spm_jobman('run', matlabbatch);
    
    Dmesh=spm_eeg_load(filename);
    
    
    %% now simulate sources on this mesh
    for s=1:Nsim,
        simpos=Dmesh.inv{1}.mesh.tess_mni.vert(simvertind(s),:); %% get location to simulate dipole on this mesh
        
        inputs=[];
        
        matlabbatch=[];
        matlabbatch{1}.spm.meeg.source.simulate.D = {filename};
        matlabbatch{1}.spm.meeg.source.simulate.val = 1;
        matlabbatch{1}.spm.meeg.source.simulate.prefix = sprintf('sim_mesh%d_source%d',simmeshind,s);
        matlabbatch{1}.spm.meeg.source.simulate.whatconditions.all = 1;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.woi = invwoi;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.isSin.foi = mean(invfoi);
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.dipmom = [10 patch_extent];
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.locs = simpos;
        matlabbatch{1}.spm.meeg.source.simulate.isSNR.setSNR = SNR;
        
        
        [a,b]=spm_jobman('run', matlabbatch);
        
        simfilename=a{1}.D{1};
        
        Dsim=spm_eeg_load(simfilename); %% simulated dataset
        
        
        %% now reconstruct onto all the meshes and look at cross val and F vals
        for meshind=1:Nmesh,
            matlabbatch=[];
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.D = {Dsim.fullfile};
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.val = 1;
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.comment = '';
            % TA % matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.mri = {'D:\MEGgroup\luzia\MQ0484.4\processed\msMQ0484-0004-00001-000208-01.img,1'};
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.mri = {'E:\MEG\s01\structural\s2014-07-11_11-51-115426-00001-00192-1_STRUCT.nii'};
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.cortex = {deblank(allmeshes(meshind,:))};
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.iskull = {''};
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.oskull = {''};
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.custom.scalp = {''};
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshres = 2;
            % TA ? % matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregdefault = {'D:\MEGgroup\luzia\MQ0484.4\processed\helmet_coreg3.mat'};
            % TA:
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).fidname = 'nas';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).specification.type = [-2.0190 110.5676 47.4948];
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).specification.type = [-83.6890 36.2986 -18.7222];
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).specification.type = [79.3990 40.6376 -17.9802];
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.useheadshape = 0;
            %
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.eeg = 'EEG BEM';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.meg = 'Single Shell';
            
            spm_jobman('run', matlabbatch);
            
            preprofilename=Dsim.fullfile;
            
            spatialmodesname=[Dsim.path filesep 'testmodes.mat'];
            [spatialmodesname,Nmodes,pctest]=spm_eeg_inv_prep_modes_xval(preprofilename, ideal_Nmodes, spatialmodesname, Nblocks,ideal_pctest);
            
            
            for methind=1:Nmeth,
                
                      
                
                %%% NOW TRY SAME WITH SPM MACHINERY
                matlabbatch=[];
                matlabbatch{1}.spm.meeg.source.invertiter.D = {preprofilename};;
                matlabbatch{1}.spm.meeg.source.invertiter.val = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.whatconditions.all = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invtype = deblank(methodnames(methind,:)); %;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.woi = invwoi;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.foi = invfoi;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.hanning = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedfile = {patchfilename}; % '<UNDEFINED>';
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedrows = 1; %'<UNDEFINED>';
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm =[patch_extent./10]; %% NB A fiddle here- need to properly quantify
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
                matlabbatch{1}.spm.meeg.source.invertiter.crossval = [pctest Nblocks];
                
                
                [a1,b1]=spm_jobman('run', matlabbatch);
                Drecon=spm_eeg_load(preprofilename);
                
                allcrosserr(simmeshind,s,meshind,methind,:)=Drecon.inv{1}.inverse.crosserr;
                allcrossF(simmeshind,s,meshind,methind,:)=Drecon.inv{1}.inverse.crossF;
                
                
            end; % for methind
            
            
        end; %% for reconstruction mesh (meshind)
        close all;
    end; % for s (sources)
end; % for simmeshind (simulatiom mesh)

% TA notes:
% loop order - sim_meshes (2), sources(8), meshes(2)(, methods(1))
% loop flow - create inversion state with mesh 1
%           - sim on specified coord on mesh 1
%           - reconstruct on mesh 1 and mesh2
%           (- etc for all coords)
%           (- repeat for mesh 2)
% allcrossF - sim_mesh(2), sim(8), mesh(2), method(1), iter(4)


methind=1;
figure(1);clf;
figure(2);clf;

for simmeshind=1:Nmesh,
    otherind=setxor(simmeshind,1:Nmesh); %% other mesh index (assuming there are just 2 meshes)
    truotherF=squeeze(allcrossF(simmeshind,:,simmeshind,methind,:)-allcrossF(simmeshind,:,otherind,methind,:)); %% reconstructed on true - reconstructed on other
    figure(1);
    subplot(Nmesh,1,simmeshind);
    %ylabel('fraction correct');
    
    
    bar(truotherF)
    ylabel('free energy');
    meanF((simmeshind-1)*Nsim+1:simmeshind*Nsim)=mean(truotherF,2)';
    
    title(sprintf('Free energy, sim mesh %d',simmeshind));
    
    %% 
    figure(2);
    truotherXval=squeeze(allcrosserr(simmeshind,:,otherind,methind,:)-allcrosserr(simmeshind,:,simmeshind,methind,:)); %%  reconstructed on other -reconstructed on true ( should be positive)
    subplot(Nmesh,1,simmeshind);
    
    bar(truotherXval)
    ylabel('X val error');
    meanXval(simmeshind)=mean(mean(truotherXval));
    title(sprintf('Cross val, sim mesh %d',simmeshind));
    correctXval(simmeshind)=length(find(truotherXval>0))/length(truotherXval(:));
end;


%% do the stats on the group
lme=[zeros(Nsim*Nmesh,1) meanF'] %% difference between true and other mesh
spm_BMS(lme,[],1)





%error('stop');


%% now set up spatial modes file with xval parameters



figure;
plot(squeeze(allcrossF(:,1,:))',squeeze(allcrosserr(:,1,:))','x'); % bat & '
hold on;
plot(squeeze(allcrossF(:,2,:))',squeeze(allcrosserr(:,2,:))','o'); % bat & '




% TA % crossval=zeros(Nmodels,Nmodels);
crossval=zeros(Nmesh,Nmesh);

% TA % Fdiff=zeros(Nmodels,Nmodels);
Fdiff=zeros(Nmesh,Nmesh);
Flmisdiff=[];
for j=1:Nmesh, % Nmodels
    for k=1:Nmesh, % Nmodels
        Fdiff(j,k)=mean(allcrossF(j,1,:)-allcrossF(k,1,:)); % bat
        
        
        errdiff=allcrosserr(j,1,:)-allcrosserr(k,1,:); % bat
        
        crossval(j,k)=length(find(errdiff<0))./pctest; %Ntest %% where j is better than k
        
        
    end; %for k
end; % for j

figure; set(gcf,'Position',[488 12 319 749])
subplot(3,1,1);
imagesc(Fdiff);colorbar;
set(gca,'XTick',[1 2],'YTick',[1 2],'XTickLabel',{'mesh 1' 'mesh 2'},'YTickLabel',{'mesh 1' 'mesh 2'})
subplot(3,1,2)
imagesc(crossval);colorbar;
set(gca,'XTick',[1 2],'YTick',[1 2],'XTickLabel',{'mesh 1' 'mesh 2'},'YTickLabel',{'mesh 1' 'mesh 2'})

dum=triu(ones(size(crossval)));
flatind=find(dum);

probF =1./(1 + exp(-Fdiff(flatind)));  %% 300ms of data

figure;
plot(Fdiff(flatind),crossval(flatind),'k.','MarkerSize',24);% 'r.',Fdiff(flatind),probF,'g.');
legend('cross val','F pval');

% figure;
% plot(Fdiff(flatind),crossval2(flatind),'r.'); %% (l % Fdiff(flatind),probF,'g.');
% legend('cross val 2','F pval');

% SAVING:
ho = findobj('Type','figure');
for k = 1:length(ho), saveas(ho(k),['sim run 3 results ' num2str(k) '.emf']), end

