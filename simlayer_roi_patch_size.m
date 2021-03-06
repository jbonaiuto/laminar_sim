function simlayer_roi_patch_size(subj_info, session_num, invfoi, SNR, varargin)
% SIMLAYER_ROI_PATCH_SIZE  Run simulations with ROI 
%   analysis, with multiple simulate/reconstruct patch size combinations
%
% Use as
%   simlayer_roi_patch_size(subjects(1), 1, [10 30], -20)
% where the first argument is the subject info structure (from create_subjects),
% the second is the session number, the third is the frequency range, and
% the fourth is the SNR (db).
% 
%   simlayer_roi_patch_size(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * surf_dir - directory containing subject surfaces
%    * mri_dir - directory containing subject MRIs
%    * dipole_moment - 10 (default) or interger - moment of simulated
%    dipole

% Parse inputs
defaults = struct('surf_dir', 'd:\pred_coding\surf', 'mri_dir', 'd:\pred_coding\mri',...
    'dipole_moment', 10);  %define default values
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

        out_path=fullfile('c:/layer_sim/ttest_results',subj_info.subj_id,...
            num2str(session_num),...
            sprintf('f%d_%d_SNR%d_dipolemoment%d_sim%d_reconstruct%d',...
                invfoi(1), invfoi(2), SNR, params.dipole_moment,...
                sim_patch_size, reconstruct_patch_size));
        simlayer_roi( subj_info, session_num, invfoi, SNR, ...
            'surf_dir', params.surf_dir, 'mri_dir', params.mri_dir,...
            'out_path', out_path, 'dipole_moment', params.dipole_moment,...
            'sim_patch_size', sim_patch_size,...
            'reconstruct_patch_size', reconstruct_patch_size);
    end
end
