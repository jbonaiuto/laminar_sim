function simlayer_patch_size(subj_info, session_num, invfoi, SNR, varargin)

% Parse inputs
defaults = struct('surf_dir', 'd:\pred_coding\surf', 'mri_dir', 'd:\pred_coding\mri',...
    'dipole_moment', 10, 'nsims', 60);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults', 'EEG');
spm_jobman('initcfg'); 

patch_sizes=[5 10];

for sp=1:length(patch_sizes)
    sim_patch_size=patch_sizes(sp);
    for rp=1:length(patch_sizes)
        reconstruct_patch_size=patch_sizes(rp);

        out_file=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d_sim%d_reconstruct%d.mat',invfoi(1),invfoi(2),SNR,params.dipole_moment,sim_patch_size,reconstruct_patch_size);

        simlayer_free_energy(subj_info, session_num, invfoi, SNR, 'surf_dir', params.surf_dir, 'mri_dir', params.mri_dir, 'out_file', out_file, 'dipole_moment', params.dipole_moment, 'sim_patch_size', sim_patch_size, 'reconstruct_patch_size', reconstruct_patch_size, 'nsims', params.nsims);

    end
end
