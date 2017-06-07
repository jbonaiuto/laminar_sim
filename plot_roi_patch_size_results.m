function plot_roi_patch_size_results(subj_info, session_num, freq, snrs,...
    varargin)
% PLOT_ROI_PATCH_SIZE_RESULTS  Plot bias and accuracy for ROI
% patch size simulations
%
% Use as
%   plot_roi_patch_size_results(subj_info, 1, [10 30], [-100 -50 -20 -5 0 5])
% where the first argument is the subject info structure (from create_subjects),
% the second is the session numner, the third is the frequency range of the
% simulated data (Hz), and the fourth is a vector of SNRs (db)
% 
%   plot_roi_patch_size_results(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces
%    * methind - 1 (default) method index: 1=EBB, 2=IID, 3=COH, 4=MSP

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10,...
    'surf_dir', 'd:\pred_coding\surf', 'methind', 1);  %define default values
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
method=methodnames{params.methind};

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

% Sim/reconstruct patch size combinations
sim_patch_sizes=[5 5 10 10];
reconstruct_patch_sizes=[5 10 5 10];

roi_correct_thresholded=zeros(length(sim_patch_sizes),length(snrs)+1,...
    params.nsims*Nmesh);
roi_correct_unthresholded=zeros(length(sim_patch_sizes),length(snrs)+1,...
    params.nsims*Nmesh);
roi_correct_significant=zeros(length(sim_patch_sizes),length(snrs)+1,...
    params.nsims*Nmesh);
roi_pial_unthresholded=zeros(length(sim_patch_sizes),length(snrs)+1,...
    params.nsims*Nmesh);
roi_pial_significant=zeros(length(sim_patch_sizes),length(snrs)+1,...
    params.nsims*Nmesh);

% Load noise (no signal) results
for idx=1:length(sim_patch_sizes)
    sim_patch_size=sim_patch_sizes(idx);
    reconstruct_patch_size=reconstruct_patch_sizes(idx);

    % Load t stats
    fname=sprintf('f%d_%d_SNR5_dipolemoment0_sim%d_reconstruct%d', ...
        freq(1), freq(2), sim_patch_size, reconstruct_patch_size);
    data_dir=fullfile('C:\layer_sim\ttest_results', subj_info.subj_id,...
        num2str(session_num), fname);            
    wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
        white_mesh, orig_pial_mesh, orig_white_mesh);
    
    roi_pial_unthresholded(idx,1,:)=wmpial_t>0;
    roi_pial_significant(idx,1,:)=abs(wmpial_t)>t_thresh;
    
    roi_correct_thresholded(idx,1,:)=[wmpial_t(1:params.nsims)<-t_thresh; wmpial_t(params.nsims+1:Nmesh*params.nsims)>t_thresh];
    roi_correct_unthresholded(idx,1,:)=[wmpial_t(1:params.nsims)<0; wmpial_t(params.nsims+1:Nmesh*params.nsims)>0];
    roi_correct_significant(idx,1,:)=abs(wmpial_t)>t_thresh;
end
    
% Load SNR results
for s=1:length(snrs)
    snr=snrs(s);
    for idx=1:length(sim_patch_sizes)
        sim_patch_size=sim_patch_sizes(idx);
        reconstruct_patch_size=reconstruct_patch_sizes(idx);

        % Load t stats
        fname=sprintf('f%d_%d_SNR%d_dipolemoment%d_sim%d_reconstruct%d',...
            freq(1), freq(2), snr, params.dipole_moment, sim_patch_size, reconstruct_patch_size);
        data_dir=fullfile('C:\layer_sim\ttest_results', subj_info.subj_id,...
            num2str(session_num), fname);  
        wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
            white_mesh, orig_pial_mesh, orig_white_mesh);
                
        roi_pial_unthresholded(idx,s+1,:)=wmpial_t>0;
        roi_pial_significant(idx,s+1,:)=abs(wmpial_t)>t_thresh;
                
        roi_correct_thresholded(idx,s+1,:)=[wmpial_t(1:params.nsims)<-t_thresh; wmpial_t(params.nsims+1:Nmesh*params.nsims)>t_thresh];
        roi_correct_unthresholded(idx,s+1,:)=[wmpial_t(1:params.nsims)<0; wmpial_t(params.nsims+1:Nmesh*params.nsims)>0];
        roi_correct_significant(idx,s+1,:)=abs(wmpial_t)>t_thresh;
    end
end

styles={'-','-','--','--'};
colors={'b','r','r','b'};
figure();
hold all;
for i=1:length(sim_patch_sizes)
    perc_correct_unthresholded_stderr=squeeze(std(roi_correct_unthresholded(i,:,:),[],3))./sqrt(Nmesh*params.nsims).*100;
    perc_correct_unthresholded=squeeze(mean(roi_correct_unthresholded(i,:,:),3)).*100.0;
    perc_correct_significant=squeeze(mean(roi_correct_significant(i,:,:),3));
    
    switch colors{i}
        case 'b'
            load('fadedblue_map');
            cm=fadedblue_map;            
        case 'r'
            load('fadedred_map');
            cm=fadedred_map;
    end
    plot_fading_line([1.25*snrs(1) snrs], perc_correct_unthresholded, ...
        perc_correct_unthresholded_stderr, perc_correct_significant, cm, styles{i});   
end
xlabel('SNR (db');
ylim([0 105]);
ylabel('% correct');

figure();
hold all;
for i=1:length(sim_patch_sizes)
    perc_pial_unthresholded_stderr=squeeze(std(roi_pial_unthresholded(i,:,:),[],3))./sqrt(Nmesh*params.nsims).*100;
    perc_pial_unthresholded=squeeze(mean(roi_pial_unthresholded(i,:,:),3)).*100.0;
    perc_pial_significant=squeeze(mean(roi_correct_significant(i,:,:),3));
    
    switch colors{i}
        case 'b'
            load('fadedblue_map');
            cm=fadedblue_map;            
        case 'r'
            load('fadedred_map');
            cm=fadedred_map;
    end
    plot_fading_line([1.25*snrs(1) snrs], perc_pial_unthresholded, ...
        perc_pial_unthresholded_stderr, perc_pial_significant, cm, styles{i});  
end
xlabel('SNR (db');
ylim([0 105]);
ylabel('% pial');

disp('Comparing classifiers');
for s=1:length(snrs)
    disp(sprintf('SNR=%df, Sim patch size 5mm, Reconstruct 5mm-Reconstruct 10mm',snrs(s)));
    both_wrong=sum(~squeeze(roi_correct_thresholded(1,s,:)) & ~squeeze(roi_correct_thresholded(2,s,:)));
    right_wrong_not_over=sum(~squeeze(roi_correct_thresholded(1,s,:)) & squeeze(roi_correct_thresholded(2,s,:)));
    over_wrong_not_right=sum(squeeze(roi_correct_thresholded(1,s,:)) & ~squeeze(roi_correct_thresholded(2,s,:)));
    both_right=sum(squeeze(roi_correct_thresholded(1,s,:)) & squeeze(roi_correct_thresholded(2,s,:)));
    McNemarextest([both_wrong,right_wrong_not_over,over_wrong_not_right,both_right],2,0.05);
    
    disp(sprintf('SNR=%df, Sim patch size 10mm, Reconstruct 10mm-Reconstruct 5mm',snrs(s)));
    both_wrong=sum(~squeeze(roi_correct_thresholded(4,s,:)) & ~squeeze(roi_correct_thresholded(3,s,:)));
    right_wrong_not_under=sum(~squeeze(roi_correct_thresholded(4,s,:)) & squeeze(roi_correct_thresholded(3,s,:)));
    under_wrong_not_right=sum(squeeze(roi_correct_thresholded(4,s,:)) & ~squeeze(roi_correct_thresholded(3,s,:)));
    both_right=sum(squeeze(roi_correct_thresholded(4,s,:)) & squeeze(roi_correct_thresholded(3,s,:)));
    McNemarextest([both_wrong,right_wrong_not_under,under_wrong_not_right,both_right],2,0.05);
end
