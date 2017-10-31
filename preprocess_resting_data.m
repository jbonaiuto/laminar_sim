function preprocess_resting_data(varargin)

% Parse inputs
defaults = struct('convert', true, 'epoch', true, 'bc', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

if params.convert
    spm_jobman('initcfg');
    matlabbatch{1}.spm.meeg.convert.dataset = {'C:\layer_sim\gb070167_sofienormative_20150506_04.ds\gb070167_sofienormative_20150506_04.meg4'};
    matlabbatch{1}.spm.meeg.convert.mode.continuous.readall = 1;
    matlabbatch{1}.spm.meeg.convert.channels{1}.all = 'all';
    matlabbatch{1}.spm.meeg.convert.outfile = 'C:\layer_sim\gb070167_sofienormative_20150506_04.ds\spmeeg_gb070167_sofienormative_20150506_04.mat';
    matlabbatch{1}.spm.meeg.convert.eventpadding = 0;
    matlabbatch{1}.spm.meeg.convert.blocksize = 3276800;
    matlabbatch{1}.spm.meeg.convert.checkboundary = 1;
    matlabbatch{1}.spm.meeg.convert.saveorigheader = 0;
    matlabbatch{1}.spm.meeg.convert.inputformat = 'ctf_old';
    matlabbatch{2}.spm.meeg.preproc.downsample.D(1) = cfg_dep('Conversion: Converted Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
    matlabbatch{2}.spm.meeg.preproc.downsample.fsample_new = 250;
    matlabbatch{2}.spm.meeg.preproc.downsample.method = 'resample';
    matlabbatch{2}.spm.meeg.preproc.downsample.prefix = 'd';
    spm_jobman('run',matlabbatch);
    
    D=spm_eeg_load('C:\layer_sim\gb070167_sofienormative_20150506_04.ds\dspmeeg_gb070167_sofienormative_20150506_04.mat');
    start_idx=min(find(abs(D(316,:)-0.0749)<0.0001));
    end_idx=start_idx+max(find(abs(D(316,start_idx:size(D,2))-0.0749)<0.0001));
    start_time=D.time(start_idx)*1000;
    end_time=D.time(end_idx)*1000;
    
    spm_jobman('initcfg');
    matlabbatch={};
    matlabbatch{1}.spm.meeg.preproc.crop.D = {'C:\layer_sim\gb070167_sofienormative_20150506_04.ds\dspmeeg_gb070167_sofienormative_20150506_04.mat'};
    matlabbatch{1}.spm.meeg.preproc.crop.timewin = [start_time end_time];
    matlabbatch{1}.spm.meeg.preproc.crop.freqwin = [-Inf Inf];
    matlabbatch{1}.spm.meeg.preproc.crop.channels{1}.all = 'all';
    matlabbatch{1}.spm.meeg.preproc.crop.prefix = 'p';
    spm_jobman('run',matlabbatch);
end

if params.epoch
    spm_jobman('initcfg');
    matlabbatch={};
    matlabbatch{1}.spm.meeg.preproc.epoch.D = {'C:\layer_sim\gb070167_sofienormative_20150506_04.ds\pdspmeeg_gb070167_sofienormative_20150506_04.mat'};
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.arbitrary.trialength = 5000;
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.arbitrary.conditionlabel = 'all';
    matlabbatch{1}.spm.meeg.preproc.epoch.bc = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.eventpadding = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.prefix = 'e';
    spm_jobman('run',matlabbatch);
end

if params.bc
    spm_jobman('initcfg');
    matlabbatch={};
    matlabbatch{1}.spm.meeg.preproc.bc.D = {'C:\layer_sim\gb070167_sofienormative_20150506_04.ds\epdspmeeg_gb070167_sofienormative_20150506_04.mat'};
    matlabbatch{1}.spm.meeg.preproc.bc.timewin = [0 5000];
    matlabbatch{1}.spm.meeg.preproc.bc.prefix = 'b';
    spm_jobman('run',matlabbatch);
end