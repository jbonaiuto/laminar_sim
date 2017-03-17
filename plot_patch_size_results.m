function plot_patch_size_results(subj_info, session_num, freq, snr, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10, 'surf_dir', 'd:\pred_coding\surf',...
    'mri_dir', 'd:\pred_coding\mri', 'sim_patch_size',0, 'reconstruct_patch_size',0);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
Nmeth=length(methodnames);
allmeshes=strvcat(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii'),...
    fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii'));
Nmesh=size(allmeshes,1);

sim_patch_sizes=[5 5 10 10];
reconstruct_patch_sizes=[5 10 5 10];

figure();
for idx=1:length(sim_patch_sizes)
    sim_patch_size=sim_patch_sizes(idx);
    reconstruct_patch_size=reconstruct_patch_sizes(idx);
    
    subplot(length(sim_patch_sizes)/2,length(reconstruct_patch_sizes)/2,idx);
    hold all;
    legend_labels={};
    for methind=1:Nmeth,       
        method=methodnames{methind};
        scores=[];
        for simmeshind=1:Nmesh,    
            data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d_sim%d_reconstruct%d.mat',freq(1),freq(2),snr,params.dipole_moment, sim_patch_size, reconstruct_patch_size));
            load(data_file);

            % F reconstructed on true - reconstructed on other
            % num simulations x number of folds
            pialWhiteF=squeeze(allcrossF(simmeshind,1:params.nsims,2,methind)-allcrossF(simmeshind,1:params.nsims,1,methind));
            scores(end+1:end+length(pialWhiteF))=pialWhiteF;     
        end        

        labels=[-1*ones(params.nsims,1); 1*ones(params.nsims,1)];
        [x,y,t,auc]=perfcurve(labels,scores,'1');
        plot(x,y);
        norm_p_auc=0;
        if min(t)<-0.01 && max(t)>0.01
            zdists=sqrt((t-0).^2);
            thresh_idx=find(zdists==min(zdists));
            p_auc=trapz(x(1:thresh_idx),y(1:thresh_idx));
            norm_p_auc=p_auc/x(thresh_idx);
            h=plot(x(thresh_idx),y(thresh_idx),'ro');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end 
        %legend_labels{end+1}=sprintf('%s, AUC=%0.2f ,norm pAUC=%0.3f', method, auc, norm_p_auc);
        legend_labels{end+1}=sprintf('%s, AUC=%0.2f', method, auc);
    end
    hold off;
    legend(legend_labels);
    title(sprintf('Simulation Size=%dmm, Reconstruction Size=%dmm', sim_patch_size, reconstruct_patch_size));
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
end
