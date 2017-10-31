function plot_abs_signal_perc_correct(subj_info, session_num, freq, white_noises,...
    varargin)
% PLOT_ABS_SIGNAL_PERC_CORRECT  Plot accuracy over white noise levels for whole brain and
% ROI analysis
%
% Use as
%   plot_abs_signal_perc_correct(subjects(1),1,[10 30],[10 20 35.5656 100 6.3246e+03 2000000])
% where the first argument is the subject info structure (from create_subjects),
% the second is the session number, the third is the frequency range of the
% simulated data (Hz), and the fourth is a vector of white noise values (fT RMS)
% 
%   plot_abs_signal_perc_correct(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 20 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 20,...
    'surf_dir', 'd:\pred_coding\surf');  %define default values
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

for i=methind:length(methodnames)
    method=methodnames{methind};
    disp(method);
    
    figure();
    hold on;
    
    perc_correct_unthresholded=zeros(1,length(white_noises));
    stderr_perc_correct_unthresholded=zeros(1,length(white_noises));
    perc_correct_significant=zeros(1,length(white_noises));
    
    disp('Whole brain');
    
    for s=1:length(white_noises)
        white_noise=white_noises(s);
        data_file=fullfile('D:\layer_sim\results_abs\',subj_info.subj_id,...
            num2str(session_num),...
            sprintf('allcrossF_f%d_%d_dipolemoment%d_noise%d.mat',freq(1),freq(2),params.dipole_moment,white_noise));
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
        disp(sprintf('white noise=%.2f fT RMS, perc=%.4f, p=%.5f', white_noise, perc_correct_unthresholded(s), pout));
    end
    plot_fading_line(log10(white_noises), perc_correct_unthresholded.*100, ...
        stderr_perc_correct_unthresholded.*100, perc_correct_significant, cm, '-');       

    disp('ROI');
    perc_correct_unthresholded=zeros(1,length(white_noises));
    stderr_perc_correct_unthresholded=zeros(1,length(white_noises));
    perc_correct_significant=zeros(1,length(white_noises));
    
    for s=1:length(white_noises)
        white_noise=white_noises(s);
        data_dir=fullfile('C:\layer_sim\ttest_results_abs', subj_info.subj_id,...
            num2str(session_num),...
            sprintf('f%d_%d_dipolemoment%d_noise%d', freq(1), freq(2), params.dipole_moment, white_noise));
        
        wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
            white_mesh, orig_pial_mesh, orig_white_mesh, 'recompute_trials', false);
        
        correct_unthresholded=[wmpial_t(1:params.nsims)<0; wmpial_t(params.nsims+1:Nmesh*params.nsims)>0];
        correct_significant=[abs(wmpial_t(1:params.nsims))>t_thresh; abs(wmpial_t(params.nsims+1:Nmesh*params.nsims))>t_thresh];
        perc_correct_unthresholded(s)=mean(correct_unthresholded);
        stderr_perc_correct_unthresholded(s)=std(correct_unthresholded)/sqrt(length(correct_unthresholded));
        perc_correct_significant(s)=mean(correct_significant);
        
        pout=myBinomTest(sum(correct_unthresholded),length(correct_unthresholded),0.5,'two');
        disp(sprintf('white noise=%.2f fT RMS, perc=%.4f, p=%.5f', white_noise, perc_correct_unthresholded(s), pout));
    end
    plot_fading_line(log10(white_noises), perc_correct_unthresholded.*100, ...
        stderr_perc_correct_unthresholded.*100, perc_correct_significant, cm, '--');           
    
    hold off;
    
    xlabel('log white noise (fT)');
    ylabel('% Correct');
    ylim([20 105]);
end

