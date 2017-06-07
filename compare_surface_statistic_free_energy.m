function compare_surface_statistic_free_energy(subj_info, session_num,...
    freq, varargin)
% COMPARE_SURFACE_STATISTIC_FREE_ENERGY  Compare the correlation
% coefficients of thickness, curvature, sulcal depth and lead field norm
% with free energy (EBB)
%
% Use as
%   compare_surface_statistic_free_energy(subj_info, 1, [10 30])
% where the first argument is the subject info structure (from create_subjects),
% the second is the session numner, and the third is the frequency range of the
% simulated data (Hz)
% 
%   compare_surface_statistic_free_energy(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * snr - -20 (default) or integer - signal to noise ratio (db)
%    * dipole_moment - 10 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces

% Parse inputs
defaults = struct('nsims', 60, 'snr', -20, 'dipole_moment', 10,...
    'surf_dir', 'd:\pred_coding\surf');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Original and downsampled white matter surface
orig_white_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_white.hires.deformed.surf.gii');

% Original and downsampled pial surface
orig_pial_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_pial.hires.deformed.surf.gii');

wm=gifti(white_mesh);

statistics={'thickness','curvature','depth','lead_field_norm'};
all_stats=[];

for i=1:length(statistics)
    statistic=statistics{i};
    
    % Compute surface statistics
    switch statistic
        case 'thickness'
            [pial_statistic, wm_statistic]=compute_thickness(pial_mesh, ...
                white_mesh, orig_pial_mesh, orig_white_mesh);
        case 'curvature'
            pial_statistic=compute_curvature(pial_mesh,...
                'curvature_type', 'mean');
            wm_statistic=compute_curvature(white_mesh,...
                'curvature_type', 'mean');
        case 'depth'
            [pial_statistic,HS]=compute_sulcal_depth(pial_mesh);
            mapping=dsearchn(HS.vertices,wm.vertices);
            wm_statistic=sqrt(sum((wm.vertices-HS.vertices(mapping,:)).^2,2));
        case 'lead_field_norm'
            D=spm_eeg_load('C:\Users\jbonai\Dropbox\meg\layer_sim\rgb_1_pial.mat');
            pial_statistic=sqrt(sum(D.inv{1}.inverse.L.^2,1))';

            D=spm_eeg_load('C:\Users\jbonai\Dropbox\meg\layer_sim\rgb_1_white.mat');
            wm_statistic=sqrt(sum(D.inv{1}.inverse.L.^2,1))';
    end
    
    % Get statistics for simulated vertices
    nverts=size(pial_statistic,1);
    rng(0);
    simvertind=randperm(nverts);
    all_stats(:,end+1)=[wm_statistic(simvertind(1:params.nsims));...
        pial_statistic(simvertind(1:params.nsims))];
end
    
Nmesh=2;
% Method = EBB
methind=1;

% Ftrue-Fother for each simulation source
allTrueOtherF=zeros(Nmesh*params.nsims,1);

% Load whole brain results
fname=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),...
    freq(2),params.snr,params.dipole_moment);
data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id,...
    num2str(session_num), fname);
load(data_file);

% Compute Ftrue - Fother for each simulation - whole brain analysis
for simmeshind=1:Nmesh,    
    % F reconstructed on true - reconstructed on other
    % num simulations x number of folds
    trueF=squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind));
    otherF=squeeze(allcrossF(simmeshind,1:params.nsims,(2-simmeshind)+1,methind));
    trueotherF=trueF-otherF;
    allTrueOtherF((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=trueotherF;
end

% Get correlation matrix for all surface stats and Ftrue-Fother
R=corr([all_stats allTrueOtherF],'type','Spearman');
% Meng's test for correlated correlation coefficients - run on abs(R)
% because we are interested in the strength of the correlation
[h,p,z]=mengz(abs(R), 5, Nmesh*params.nsims);
disp(sprintf('Chi-sq=%.3f, p=%.5f', z, p));

% Pairwise comparisons
for i=1:length(statistics)
    stat_one=statistics{i};
    for j=i+1:length(statistics)
        stat_two=statistics{j};
        lambda=zeros(1,4);
        lambda(i)=1.0;
        lambda(j)=-1.0;
        [h,p,z]=mengz(abs(R),5,Nmesh*params.nsims,lambda);
        disp(sprintf('Comparison: %s>%s, Z=%.2f, p=%.4f', stat_one, stat_two, z, p));
        
        lambda(i)=-1.0;
        lambda(j)=1.0;
        [h,p,z]=mengz(abs(R),5,Nmesh*params.nsims,lambda);
        disp(sprintf('Comparison: %s<%s, Z=%.2f, p=%.4f', stat_one, stat_two, z, p));
    end
end
