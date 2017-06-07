function compare_classifier_accuracy(subj_info, session_num, freq, snr,...
    varargin)
% COMPARE_CLASSIFIER_ACCURACY  Compare the accuracy of ROI and whole brain
%   analyses using multiple inversion algorithms
%
% Use as
%   compare_classifier_accuracy(subj_info, 1, [10 30], -20)
% where the first argument is the subject info structure (from create_subjects),
% the second is the session numner, the third is the frequency range of the
% simulated data (Hz), and the fourth is the SNR (db)
% 
%   compare_classifier_accuracy(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10, ...
    'surf_dir', 'd:\pred_coding\surf');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Critical t-statistic
dof=514;
alpha=1.0-0.05/2;
t_thresh=tinv(alpha, dof);

methodnames={'EBB','IID','COH','MSP'};
Nmeth=length(methodnames);

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

simmeshes={white_mesh,pial_mesh};
Nmesh=length(simmeshes);

% Correct classifications for each method and simulation - whole brain (F)
% and ROI (t)
f_correct=zeros(Nmeth,Nmesh*params.nsims);
t_correct=zeros(Nmeth,Nmesh*params.nsims);

% Load whole brain results
fname=sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date);
data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, ...
    num2str(session_num), fname);
load(data_file);

% Compute Ftrue - Fother for each simulation - whole brain analysis
for methind=1:Nmeth    
    for simmeshind=1:Nmesh,    
        % F value for simulated surface
        trueF=squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind));
        % F value for other surface
        otherF=squeeze(allcrossF(simmeshind,1:params.nsims,(2-simmeshind)+1,methind));
        trueOtherF=trueF-otherF;
        % Correct where Ftrue-Fother>3
        for x=1:params.nsims
            f_correct(methind,(simmeshind-1)*params.nsims+x)=trueOtherF(x)>3;
        end
    end
end

% Get t-stat for each simulation - ROI analysis
for methind=1:Nmeth
    method=methodnames{methind};
    fname=sprintf('f%d_%d_SNR%d_dipolemoment%d', freq(1), freq(2), snr,...
        params.dipole_moment);
    data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id,...
        num2str(session_num), fname);
        
    % Get wm and pial simulation t-statistics
    wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
        white_mesh, orig_pial_mesh, orig_white_mesh,...
        'recompute_trials', false);

    % Correct when white matter simulation and t<-t_thresh or pial
    % simulation and t>t_thresh
    for simmeshind=1:Nmesh,    
        for s=1:params.nsims
            tstat=wmpial_t((simmeshind-1)*params.nsims+s);
            t_correct(methind,(simmeshind-1)*params.nsims+s)=(simmeshind==1 && tstat<-t_thresh) || (simmeshind==2 && tstat>t_thresh);            
        end
    end  
end

disp('Comparing classifiers');

% McNemar test - EBB (whole brain) - MSP (whole brain)
disp('EBB (whole brain) - MSP (whole brain)');
both_wrong=sum(~f_correct(1,:) & ~f_correct(4,:));
ebb_wrong_not_msp=sum(~f_correct(1,:) & f_correct(4,:));
msp_wrong_not_ebb=sum(f_correct(1,:) & ~f_correct(4,:));
both_right=sum(f_correct(1,:) & f_correct(4,:));
McNemarextest([both_wrong,ebb_wrong_not_msp,msp_wrong_not_ebb,both_right],2,0.05);

% McNemar test - EBB (ROI) - MSP (ROI)
disp('ROI, EBB-MSP');
both_wrong=sum(~t_correct(1,:) & ~t_correct(4,:));
ebb_wrong_not_msp=sum(~t_correct(1,:) & t_correct(4,:));
msp_wrong_not_ebb=sum(t_correct(1,:) & ~t_correct(4,:));
both_right=sum(t_correct(1,:) & t_correct(4,:));
McNemarextest([both_wrong,ebb_wrong_not_msp,msp_wrong_not_ebb,both_right],2,0.05);

% McNemar test - Whole brain (EBB) - ROI (EBB)
disp('EBB, Whole brain-ROI:');
both_wrong=sum(~f_correct(1,:) & ~t_correct(1,:));
wb_wrong_not_roi=sum(~f_correct(1,:) & t_correct(1,:));
roi_wrong_not_wb=sum(f_correct(1,:) & ~t_correct(1,:));
both_right=sum(f_correct(1,:) & t_correct(1,:));
McNemarextest([both_wrong,wb_wrong_not_roi,roi_wrong_not_wb,both_right],2,0.05);

% McNemar test - Whole brain (MSP) - ROI (MSP)
disp('MSP, Whole brain-ROI:');
both_wrong=sum(~f_correct(4,:) & ~t_correct(4,:));
wb_wrong_not_roi=sum(~f_correct(4,:) & t_correct(4,:));
roi_wrong_not_wb=sum(f_correct(4,:) & ~t_correct(4,:));
both_right=sum(f_correct(4,:) & t_correct(4,:));
McNemarextest([both_wrong,wb_wrong_not_roi,roi_wrong_not_wb,both_right],2,0.05);

