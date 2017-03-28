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

labels=[-1*ones(params.nsims,1); 1*ones(params.nsims,1)];
f_scores=zeros(length(methodnames),length(sim_patch_sizes),params.nsims*Nmesh);

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
        f_scores(methind,idx,:)=scores;
        
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

labels(find(labels<0))=0;
[pvalue Wft Wtf]=pauc(squeeze(f_scores(4,2,:)),squeeze(f_scores(4,3,:)),labels);
disp(sprintf('MSP, S5R10-S10R5, W_5-10=%.4f, W_10-5=%.4f, p=%.5f', Wft, Wtf, pvalue));

[pvalue Wft Wtf]=pauc(squeeze(f_scores(1,2,:)),squeeze(f_scores(1,3,:)),labels);
disp(sprintf('EBB, S5R10-S10R5, W_5-10=%.4f, W_10-5=%.4f, p=%.5f', Wft, Wtf, pvalue));

[pvalue Webb Wmsp]=pauc(squeeze(f_scores(1,2,:)),squeeze(f_scores(4,2,:)),labels);
disp(sprintf('S5R10, EBB-MSP, W_EBB=%.4f, W_MSP=%.4f, p=%.5f', Webb, Wmsp, pvalue));

[pvalue Webb Wmsp]=pauc(squeeze(f_scores(1,3,:)),squeeze(f_scores(4,3,:)),labels);
disp(sprintf('S10R5, EBB-MSP, W_EBB=%.4f, W_MSP=%.4f, p=%.5f', Webb, Wmsp, pvalue));
