function plot_snr_perc_pial(subj_info, session_num, freq, snrs, varargin)
% PLOT_SNR_PERC_PIAL  Plot pial bias over SNR levels for whole brain and
% ROI analysis
%
% Use as
%   plot_snr_perc_pial(subj_info, 1, [10 30], [-100 -50 -20 -5 0 5])
% where the first argument is the subject info structure (from create_subjects),
% the second is the session numner, the third is the frequency range of the
% simulated data (Hz), and the fourth is a vector of SNRs (db)
% 
%   plot_snr_perc_pial(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10,...
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
    
    perc_pial_unthresholded=zeros(1,length(snrs)+1);
    stderr_perc_pial_unthresholded=zeros(1,length(snrs)+1);
    perc_pial_significant=zeros(1,length(snrs)+1);
    
    disp('whole brain');
    
    % Noise (-inf db)
    data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id,...
        num2str(session_num),...
        sprintf('allcrossF_f%d_%d_SNR5_dipolemoment0.mat',freq(1),freq(2)));
    load(data_file);

    pial_unthresholded=zeros(1,Nmesh*params.nsims);
    pial_significant=zeros(1,Nmesh*params.nsims);
    for simmeshind=1:Nmesh,    
        pialwmF=squeeze(allcrossF(simmeshind,1:params.nsims,2,methind)-allcrossF(simmeshind,1:params.nsims,1,methind));
        pial_unthresholded((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=pialwmF>0;
        pial_significant((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=abs(pialwmF)>3;
    end
    perc_pial_unthresholded(1)=mean(pial_unthresholded);
    stderr_perc_pial_unthresholded(1)=std(pial_unthresholded)/sqrt(length(pial_unthresholded));
    perc_pial_significant(1)=mean(pial_significant);
        
    pout=myBinomTest(sum(pial_unthresholded),length(pial_unthresholded),0.5,'two');
    disp(sprintf('SNR=-inf, pial=%.2f, p=%.5f', perc_pial_unthresholded(1)*100.0, pout));
    
    for s=1:length(snrs)
        snr=snrs(s);
        data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id,...
            num2str(session_num),...
            sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),snr, params.dipole_moment));
        load(data_file);

        pial_unthresholded=zeros(1,Nmesh*params.nsims);
        pial_significant=zeros(1,Nmesh*params.nsims);
        for simmeshind=1:Nmesh,    
            pialwmF=squeeze(allcrossF(simmeshind,1:params.nsims,2,methind)-allcrossF(simmeshind,1:params.nsims,1,methind));
            pial_unthresholded((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=pialwmF>0;
            pial_significant((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=abs(pialwmF)>3;            
        end
        perc_pial_unthresholded(s+1)=mean(pial_unthresholded);
        stderr_perc_pial_unthresholded(s+1)=std(pial_unthresholded)/sqrt(length(pial_unthresholded));
        perc_pial_significant(s+1)=mean(pial_significant);
        
        pout=myBinomTest(sum(pial_unthresholded),length(pial_unthresholded),0.5,'two');
        disp(sprintf('SNR=%d dB, pial=%.2f, p=%.5f', snr, perc_pial_unthresholded(s+1)*100.0, pout));
    end
    plot_fading_line([1.25*snrs(1) snrs], perc_pial_unthresholded.*100, ...
        stderr_perc_pial_unthresholded.*100, perc_pial_significant, cm, '-');
    
    disp('ROI');
    perc_pial_unthresholded=zeros(1,length(snrs)+1);
    stderr_perc_pial_unthresholded=zeros(1,length(snrs)+1);
    perc_pial_significant=zeros(1,length(snrs)+1);
    
    % Noise (-inf db)
    data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id,...
        num2str(session_num), ...
        sprintf('f%d_%d_SNR5_dipolemoment0', freq(1), freq(2)));
        
    wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
        white_mesh, orig_pial_mesh, orig_white_mesh);

    pial_unthresholded=wmpial_t>0;
    pial_significant=abs(wmpial_t)>t_thresh;
    perc_pial_unthresholded(1)=mean(pial_unthresholded);
    perc_pial_significant(1)=mean(pial_significant);
    stderr_perc_pial_unthresholded(1)=std(pial_unthresholded)/sqrt(length(pial_unthresholded));
    
    pout=myBinomTest(sum(pial_unthresholded),length(pial_unthresholded),0.5,'two');
    disp(sprintf('SNR=-inf, pial=%.2f, p=%.5f', perc_pial_unthresholded(1)*100.0, pout));
        
    for s=1:length(snrs)
        snr=snrs(s);
        data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id,...
            num2str(session_num),...
            sprintf('f%d_%d_SNR%d_dipolemoment%d', freq(1), freq(2), snr, params.dipole_moment));
        
        wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
            white_mesh, orig_pial_mesh, orig_white_mesh);
        
        pial_unthresholded=wmpial_t>0;      
        pial_significant=abs(wmpial_t)>t_thresh;
        perc_pial_unthresholded(s+1)=mean(pial_unthresholded);
        stderr_perc_pial_unthresholded(s+1)=std(pial_unthresholded)/sqrt(length(pial_unthresholded));
        perc_pial_significant(s+1)=mean(pial_significant);
        
        pout=myBinomTest(sum(pial_unthresholded),length(pial_unthresholded),0.5,'two');
        disp(sprintf('SNR=%d dB, pial=%.2f, p=%.5f', snr, perc_pial_unthresholded(s+1)*100.0, pout));
    end
    plot_fading_line([1.25*snrs(1) snrs], perc_pial_unthresholded.*100, ...
        stderr_perc_pial_unthresholded.*100, perc_pial_significant, cm, '--');
    
    hold off;
    
    xlabel('SNR');
    ylabel('% Pial');
    ylim([0 105]);
end

