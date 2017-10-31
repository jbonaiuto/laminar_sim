function plot_snr_perc_correct_ve(subj_info, session_num, freq, snrs,...
    varargin)
% PLOT_SNR_PERC_CORRECT_VE  Plot accuracy over SNR levels for whole brain analysis
% using variance explained rather than free energy
%
% Use as
%   plot_snr_perc_correct_ve(subj_info, 1, [10 30], [-100 -50 -20 -5 0 5])
% where the first argument is the subject info structure (from create_subjects),
% the second is the session numner, the third is the frequency range of the
% simulated data (Hz), and the fourth is a vector of SNRs (db)
% 
%   plot_snr_perc_correct_ve(...,'param','value','param','value'...) allows
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

methodnames={'EBB','IID','COH','MSP'};

% Downsampled white matter surface
white_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_white.hires.deformed.surf.gii');

% Downsampled pial surface
pial_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_pial.hires.deformed.surf.gii');

simmeshes={white_mesh,pial_mesh};
Nmesh=length(simmeshes);

for i=methind:length(methodnames)
    method=methodnames{methind};
    disp(method);
    
    figure();
    title(method);
    hold all;
    
    perc_correct_unthresholded=zeros(1,length(snrs)+1);
    stderr_perc_correct_unthresholded=zeros(1,length(snrs)+1);
        
    % Noise (-inf db)
    data_file=fullfile('C:\layer_sim\results\',subj_info.subj_id,...
        num2str(session_num),sprintf('allcrossF_f%d_%d_SNR5_dipolemoment0.mat',freq(1),freq(2)));
    load(data_file);

    correct_unthresholded=zeros(1,Nmesh*params.nsims);
    for simmeshind=1:Nmesh,    
        truotherVE=squeeze(allcrossVE(simmeshind,1:params.nsims,simmeshind,methind)-allcrossVE(simmeshind,1:params.nsims,(2-simmeshind)+1,methind));
        correct_unthresholded((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=truotherVE>0;
    end
    perc_correct_unthresholded(1)=mean(correct_unthresholded);
    stderr_perc_correct_unthresholded(1)=std(correct_unthresholded)/sqrt(length(correct_unthresholded));

    pout=myBinomTest(sum(correct_unthresholded),length(correct_unthresholded),0.5,'two');
    disp(sprintf('SNR=-inf, accuracy=%.2f, p=%.5f', perc_correct_unthresholded(1)*100.0, pout));

    for s=1:length(snrs)
        snr=snrs(s);
        data_file=fullfile('C:\layer_sim\results\',subj_info.subj_id,...
            num2str(session_num),sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),snr, params.dipole_moment));
        load(data_file);

        correct_unthresholded=zeros(1,Nmesh*params.nsims);
        for simmeshind=1:Nmesh,    
            truotherVE=squeeze(allcrossVE(simmeshind,1:params.nsims,simmeshind,methind)-allcrossVE(simmeshind,1:params.nsims,(2-simmeshind)+1,methind));
            correct_unthresholded((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=truotherVE>0;
        end
        perc_correct_unthresholded(s+1)=mean(correct_unthresholded);
        stderr_perc_correct_unthresholded(s+1)=std(correct_unthresholded)/sqrt(length(correct_unthresholded));

        pout=myBinomTest(sum(correct_unthresholded),length(correct_unthresholded),0.5,'two');
        disp(sprintf('SNR=%.2f dB, accuracy=%.2f, p=%.5f', snr, perc_correct_unthresholded(s+1)*100.0, pout));
    end
    errorbar([1.25*snrs(1) snrs], perc_correct_unthresholded.*100, stderr_perc_correct_unthresholded.*100,...
        'MarkerSize',6);
    
    hold off;
    xlabel('SNR');
    ylabel('% Correct');
    ylim([20 105]);
end

