function plot_real_noise_perc_correct(subj_info, session_num, freq, dipole_moments,...
    varargin)
% PLOT_REAL_NOISE_PERC_CORRECT  Plot accuracy over dipole moment levels for whole brain and
% ROI analysis where resting state data used as noise
%
% Use as
%   plot_real_noise_perc_correct(subj_info, 1, [10 30], [20 30 40 50 60])
% where the first argument is the subject info structure (from create_subjects),
% the second is the session numner, the third is the frequency range of the
% simulated data (Hz), and the fourth is a vector of dipole moments (nAm)
% 
%   plot_real_noise_perc_correct(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * surf_dir - directory containing subject surfaces

% Parse inputs
defaults = struct('nsims', 60, 'surf_dir', 'd:\pred_coding\surf');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% critical t value
dof=514;
alpha=1.0-0.05/2;
t_thresh=tinv(alpha, dof);

methodnames={'EBB','IID','COH','MSP'};

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

load('fadedblue_map');
cm=fadedblue_map;

for methind=1:length(methodnames)
    method=methodnames{methind};
    disp(method);
    
    figure();
    hold on;
    
    disp('Whole brain');
    perc_correct_unthresholded=zeros(1,length(dipole_moments));
    stderr_perc_correct_unthresholded=zeros(1,length(dipole_moments));
    perc_correct_significant=zeros(1,length(dipole_moments));
    
    for s=1:length(dipole_moments)
        dipole_moment=dipole_moments(s);
        data_file=fullfile('D:\layer_sim\results_real_noise\',subj_info.subj_id,...
            num2str(session_num),...
            sprintf('allcrossF_f%d_%d_dipolemoment%d.mat',freq(1),freq(2), dipole_moment));
        load(data_file);

        correct_unthresholded=zeros(1,Nmesh*params.nsims);
        correct_significant=zeros(1,Nmesh*params.nsims);
        for simmeshind=1:Nmesh,    
            truotherF=squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind)-allcrossF(simmeshind,1:params.nsims,(2-simmeshind)+1,methind));
            correct_unthresholded((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=truotherF>0;
            correct_significant((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=abs(truotherF)>3;
        end
        perc_correct_unthresholded(s)=mean(correct_unthresholded);
        stderr_perc_correct_unthresholded(s)=std(correct_unthresholded)/sqrt(length(correct_unthresholded));
        perc_correct_significant(s)=mean(correct_significant);
        
        pout=myBinomTest(sum(correct_unthresholded),length(correct_unthresholded),0.5,'two');
        disp(sprintf('dipole moment=%.2f nAm, perc=%.4f, p=%.5f', dipole_moment, perc_correct_unthresholded(s), pout));
    end
    plot_fading_line(dipole_moments, perc_correct_unthresholded.*100, ...
        stderr_perc_correct_unthresholded.*100, perc_correct_significant, cm, '-');       
    
    disp('ROI');
    perc_correct_unthresholded=zeros(1,length(dipole_moments));
    stderr_perc_correct_unthresholded=zeros(1,length(dipole_moments));
    perc_correct_significant=zeros(1,length(dipole_moments));
    
    for s=1:length(dipole_moments)
        dipole_moment=dipole_moments(s);
        data_dir=fullfile('D:\layer_sim\ttest_results_real_noise', subj_info.subj_id,...
            num2str(session_num),...
            sprintf('f%d_%d_dipolemoment%d', freq(1), freq(2), dipole_moment));
        
        wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
            white_mesh, orig_pial_mesh, orig_white_mesh, 'recompute',false, 'recompute_trials', false);
        
        correct_unthresholded=[wmpial_t(1:params.nsims)<0; wmpial_t(params.nsims+1:Nmesh*params.nsims)>0];
        correct_significant=[abs(wmpial_t(1:params.nsims))>t_thresh; abs(wmpial_t(params.nsims+1:Nmesh*params.nsims))>t_thresh];
        perc_correct_unthresholded(s)=mean(correct_unthresholded);
        stderr_perc_correct_unthresholded(s)=std(correct_unthresholded)/sqrt(length(correct_unthresholded));
        perc_correct_significant(s)=mean(correct_significant);
        
        pout=myBinomTest(sum(correct_unthresholded),length(correct_unthresholded),0.5,'two');
        disp(sprintf('dipole moment=%.2f nAm, perc=%.4f, p=%.5f', dipole_moment, perc_correct_unthresholded(s), pout));
    end
    plot_fading_line(dipole_moments, perc_correct_unthresholded.*100, ...
        stderr_perc_correct_unthresholded.*100, perc_correct_significant, cm, '--');           
    
    hold off;
    
    xlabel('Dipole moment (nAm)');
    ylabel('% Correct');
    ylim([20 105]);
end

