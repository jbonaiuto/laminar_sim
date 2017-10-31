function plot_abs_signal_perc_pial(subj_info, session_num, freq, white_noises, varargin)
% PLOT_ABS_SIGNAL_PERC_PIAL  Plot pial bias over white noise levels for whole brain and
% ROI analysis
%
% Use as
%   plot_abs_signal_perc_pial(subj_info, 1, [10 30], [10 20 35.5656 100 6.3246e+03 2000000])
% where the first argument is the subject info structure (from create_subjects),
% the second is the session numner, the third is the frequency range of the
% simulated data (Hz), and the fourth is a vector of white noise levels (fT RMS)
% 
%   plot_abs_signal_perc_pial(...,'param','value','param','value'...) allows
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

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now

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
    
    perc_pial_unthresholded=zeros(1,length(white_noises));
    stderr_perc_pial_unthresholded=zeros(1,length(white_noises));
    perc_pial_significant=zeros(1,length(white_noises));
    
    disp('Whole brain');
    
    for s=1:length(white_noises)
        white_noise=white_noises(s);
        data_file=fullfile('D:\layer_sim\results_abs\',subj_info.subj_id,...
            num2str(session_num),...
            sprintf('allcrossF_f%d_%d_dipolemoment%d_noise%d.mat',freq(1),freq(2),params.dipole_moment,white_noise));
        load(data_file);

        pial_unthresholded=zeros(1,Nmesh*params.nsims);
        pial_significant=zeros(1,Nmesh*params.nsims);
        for simmeshind=1:Nmesh,    
            pialwmF=squeeze(allcrossF(simmeshind,1:params.nsims,2,methind)-allcrossF(simmeshind,1:params.nsims,1,methind));
            pial_unthresholded((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=pialwmF>0;
            pial_significant((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=abs(pialwmF)>3;            
        end
        perc_pial_unthresholded(s)=mean(pial_unthresholded);
        stderr_perc_pial_unthresholded(s)=std(pial_unthresholded)/sqrt(length(pial_unthresholded));
        perc_pial_significant(s)=mean(pial_significant);
        
        pout=myBinomTest(sum(pial_unthresholded),length(pial_unthresholded),0.5,'two');
        disp(sprintf('white noise=%.2f fT RMS, perc=%.4f, p=%.5f', white_noise, perc_pial_unthresholded(s), pout));
    end
    plot_fading_line(log10(white_noises), perc_pial_unthresholded.*100, ...
        stderr_perc_pial_unthresholded.*100, perc_pial_significant, cm, '-');
    
    disp('ROI');
    perc_pial_unthresholded=zeros(1,length(white_noises));
    stderr_perc_pial_unthresholded=zeros(1,length(white_noises));
    perc_pial_significant=zeros(1,length(white_noises));
    
    for s=1:length(white_noises)
        white_noise=white_noises(s);
        data_dir=fullfile('c:\layer_sim\ttest_results_abs', subj_info.subj_id,...
            num2str(session_num),...
            sprintf('f%d_%d_dipolemoment%d_noise%d', freq(1), freq(2), params.dipole_moment, white_noise));
        
        wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
            white_mesh, orig_pial_mesh, orig_white_mesh);
        
        pial_unthresholded=wmpial_t>0;      
        pial_significant=abs(wmpial_t)>t_thresh;
        perc_pial_unthresholded(s)=mean(pial_unthresholded);
        stderr_perc_pial_unthresholded(s)=std(pial_unthresholded)/sqrt(length(pial_unthresholded));
        perc_pial_significant(s)=mean(pial_significant);
        
        pout=myBinomTest(sum(pial_unthresholded),length(pial_unthresholded),0.5,'two');
        disp(sprintf('white noise=%.2f fT RMS, perc=%.4f, p=%.5f', white_noise, perc_pial_unthresholded(s), pout));
    end
    plot_fading_line(log10(white_noises), perc_pial_unthresholded.*100, ...
        stderr_perc_pial_unthresholded.*100, perc_pial_significant, cm, '--');
    
    hold off;
    
    xlabel('White noise (fT)');
    ylabel('% Pial');
    ylim([0 105]);
end

