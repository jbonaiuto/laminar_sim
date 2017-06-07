function plot_free_energy_patch_size_results(subj_info, session_num, freq,...
    snrs, varargin)
% PLOT_FREE_ENERGY_PATCH_SIZE_RESULTS  Plot bias and accuracy for whole
% brain patch size simulations
%
% Use as
%   plot_free_energy_patch_size_results(subj_info, 1, [10 30], [-100 -50 -20 -5 0 5])
% where the first argument is the subject info structure (from create_subjects),
% the second is the session numner, the third is the frequency range of the
% simulated data (Hz), and the fourth is a vector of SNRs (db)
% 
%   compare_surface_statistic_free_energy(...,'param','value','param','value'...) allows
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

methodnames={'EBB','IID','COH','MSP'};
Nmeth=length(methodnames);

white_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_white.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_pial.hires.deformed.surf.gii');

allmeshes=strvcat(white_mesh, pial_mesh);
Nmesh=size(allmeshes,1);

% Sim/reconstruct patch size combinations
sim_patch_sizes=[5 5 10 10];
reconstruct_patch_sizes=[5 10 5 10];

f_correct_thresholded=zeros(length(methodnames),length(sim_patch_sizes),...
    length(snrs)+1,params.nsims*Nmesh);
f_correct_unthresholded=zeros(length(methodnames),length(sim_patch_sizes),...
    length(snrs)+1,params.nsims*Nmesh);
f_correct_significant=zeros(length(methodnames),length(sim_patch_sizes),...
    length(snrs)+1,params.nsims*Nmesh);
f_pial_unthresholded=zeros(length(methodnames),length(sim_patch_sizes),...
    length(snrs)+1,params.nsims*Nmesh);
f_pial_significant=zeros(length(methodnames),length(sim_patch_sizes),...
    length(snrs)+1,params.nsims*Nmesh);

% Load noise (no signal) results
for idx=1:length(sim_patch_sizes)
    sim_patch_size=sim_patch_sizes(idx);
    reconstruct_patch_size=reconstruct_patch_sizes(idx);

    % Load whole brain results - no signal (noise
    fname=sprintf('allcrossF_f%d_%d_SNR5_dipolemoment0_sim%d_reconstruct%d.mat',...
        freq(1),freq(2),sim_patch_size, reconstruct_patch_size);
    data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id,...
        num2str(session_num), fname);
    load(data_file);

    for methind=1:Nmeth,       
        for simmeshind=1:Nmesh,                
            % F reconstructed on pial - reconstructed on white matter
            pialF=squeeze(allcrossF(simmeshind,1:params.nsims,2,methind));
            whiteF=squeeze(allcrossF(simmeshind,1:params.nsims,1,methind));
            pialWhiteF=pialF-whiteF;
            f_pial_unthresholded(methind,idx,1,...
                (simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=pialWhiteF>0;  
            f_pial_significant(methind,idx,1,...
                (simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=abs(pialWhiteF)>3;  

            % F reconstructed on true - reconstructed on other
            trueF=squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind));
            otherF=squeeze(allcrossF(simmeshind,1:params.nsims,2-simmeshind+1,methind));
            trueOtherF=trueF-otherF;
            f_correct_thresholded(methind,idx,1,...
                (simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=trueOtherF>3;
            f_correct_unthresholded(methind,idx,1,...
                (simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=trueOtherF>0;
            f_correct_significant(methind,idx,1,...
                (simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=abs(trueOtherF)>3;
        end        
    end
end
    
% Load results for each SNR level
for s=1:length(snrs)
    snr=snrs(s);
    for idx=1:length(sim_patch_sizes)
        sim_patch_size=sim_patch_sizes(idx);
        reconstruct_patch_size=reconstruct_patch_sizes(idx);

        % Load whole brain results
        fname=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d_sim%d_reconstruct%d.mat',...
            freq(1),freq(2),snr,params.dipole_moment, sim_patch_size,...
            reconstruct_patch_size);
        data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id,...
            num2str(session_num), fname);
        load(data_file);

        for methind=1:Nmeth,       
            for simmeshind=1:Nmesh,                   
                % F reconstructed on pial - reconstructed on white matter
                pialF=squeeze(allcrossF(simmeshind,1:params.nsims,2,methind));
                whiteF=squeeze(allcrossF(simmeshind,1:params.nsims,1,methind));
                pialWhiteF=pialF-whiteF;
                f_pial_unthresholded(methind,idx,s+1,...
                    (simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=pialWhiteF>0; 
                f_pial_significant(methind,idx,s+1,...
                    (simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=abs(pialWhiteF)>3;  

                % F reconstructed on true - reconstructed on other
                trueF=squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind));
                otherF=squeeze(allcrossF(simmeshind,1:params.nsims,2-simmeshind+1,methind));
                trueOtherF=trueF-otherF;
                f_correct_thresholded(methind,idx,s+1,...
                    (simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=trueOtherF>3;
                f_correct_unthresholded(methind,idx,s+1,...
                    (simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=trueOtherF>0;
                f_correct_significant(methind,idx,s+1,...
                    (simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=abs(trueOtherF)>3;
            end        
        end
    end
end

styles={'-','-','--','--'};
colors={'b','r','r','b'};
f=figure();
hold all;
for i=1:length(sim_patch_sizes)
    perc_correct_unthresholded_stderr=squeeze(std(f_correct_unthresholded(params.methind,i,:,:),[],4))./sqrt(Nmesh*params.nsims).*100;
    perc_correct_unthresholded=squeeze(mean(f_correct_unthresholded(params.methind,i,:,:),4)).*100.0;
    perc_correct_significant=squeeze(mean(f_correct_significant(params.methind,i,:,:),4));
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
ylim([40 105]);
ylabel('% correct');

f=figure();
hold all;
for i=1:length(sim_patch_sizes)
    perc_pial_unthresholded_stderr=squeeze(std(f_pial_unthresholded(params.methind,i,:,:),[],4))./sqrt(Nmesh*params.nsims).*100;
    perc_pial_unthresholded=squeeze(mean(f_pial_unthresholded(params.methind,i,:,:),4)).*100.0;
    perc_pial_significant=squeeze(mean(f_pial_significant(params.methind,i,:,:),4));
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
for s=2:length(snrs)+1
    disp(sprintf('SNR=%df, Sim patch size 5mm, Reconstruct 5mm-Reconstruct 10mm',snrs(s-1)));
    both_wrong=sum(~squeeze(f_correct_thresholded(params.methind,1,s,:)) & ~squeeze(f_correct_thresholded(params.methind,2,s,:)));
    right_wrong_not_over=sum(~squeeze(f_correct_thresholded(params.methind,1,s,:)) & squeeze(f_correct_thresholded(params.methind,2,s,:)));
    over_wrong_not_right=sum(squeeze(f_correct_thresholded(params.methind,1,s,:)) & ~squeeze(f_correct_thresholded(params.methind,2,s,:)));
    both_right=sum(squeeze(f_correct_thresholded(params.methind,1,s,:)) & squeeze(f_correct_thresholded(params.methind,2,s,:)));
    McNemarextest([both_wrong,right_wrong_not_over,over_wrong_not_right,both_right],2,0.05);
    
    disp(sprintf('SNR=%df, Sim patch size 10mm, Reconstruct 10mm-Reconstruct 5mm',snrs(s-1)));
    both_wrong=sum(~squeeze(f_correct_thresholded(params.methind,4,s,:)) & ~squeeze(f_correct_thresholded(params.methind,3,s,:)));
    right_wrong_not_under=sum(~squeeze(f_correct_thresholded(params.methind,4,s,:)) & squeeze(f_correct_thresholded(params.methind,3,s,:)));
    under_wrong_not_right=sum(squeeze(f_correct_thresholded(params.methind,4,s,:)) & ~squeeze(f_correct_thresholded(params.methind,3,s,:)));
    both_right=sum(squeeze(f_correct_thresholded(params.methind,4,s,:)) & squeeze(f_correct_thresholded(params.methind,3,s,:)));
    McNemarextest([both_wrong,right_wrong_not_under,under_wrong_not_right,both_right],2,0.05);
end
