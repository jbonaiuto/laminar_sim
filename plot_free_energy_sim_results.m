function plot_free_energy_sim_results(subj_info, session_num, freq, snr, varargin)

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
allPialWhiteF=zeros(Nmeth,Nmesh,params.nsims);
figure();
mesh_idx=[3 4 7 8;1 2 5 6];
for methind=1:Nmeth,       
    for simmeshind=1:Nmesh,    
        subplot(4,2,mesh_idx(simmeshind,methind));
        [path,file,ext]=fileparts(deblank(allmeshes(simmeshind,:)));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        simmeshname=y{2};
        
        if params.sim_patch_size==0 && params.reconstruct_patch_size==0
            data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),snr,params.dipole_moment));
        else
            data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d_sim%d_reconstruct%d.mat',freq(1),freq(2),snr,params.dipole_moment, params.sim_patch_size, params.reconstruct_patch_size));
        end
        load(data_file);

        % F reconstructed on true - reconstructed on other
        % num simulations x number of folds
        pialWhiteF=squeeze(allcrossF(simmeshind,1:params.nsims,2,methind)-allcrossF(simmeshind,1:params.nsims,1,methind));
        allPialWhiteF(methind,simmeshind,:)=pialWhiteF;                        
        
        %bar(squeeze(allPialWhiteF(methind,simmeshind,:,:))')
        hold on
        for simind=1:params.nsims
            color='b';
            if pialWhiteF(simind)<0
                color='r';
            end
            bar(simind,pialWhiteF(simind),color);
        end
        plot([0 params.nsims+1],[3 3],'r--');
        plot([0 params.nsims+1],[-3 -3],'r--');
        xlim([0 params.nsims+1]);
        xlabel('Simulation')
        ylabel('Free energy diff (pial-white)');
        title(sprintf('Free energy, %s, %s',methodnames{methind},simmeshname));        
    end

end

figure();
hold all;
legend_labels={};

for methind=1:Nmeth,           
    method=methodnames{methind};

    labels=[-1*ones(params.nsims,1); 1*ones(params.nsims,1)];
    scores=zeros(params.nsims*Nmesh,1);    
    for simmeshind=1:Nmesh,    
        scores((simmeshind-1)*params.nsims+1:simmeshind*params.nsims,:)=squeeze(allPialWhiteF(methind,simmeshind,:));
    end
    [x,y,t,auc]=compute_roc(scores,labels,params.nsims);
    plot(x,y);
    norm_p_auc=0;
    if min(t)<0 && max(t)>0
        thresh_idx=find(t==0);
        p_auc=trapz(x(1:thresh_idx),y(1:thresh_idx));
        norm_p_auc=p_auc/x(thresh_idx);
        h=plot(x(thresh_idx),y(thresh_idx),'o','MarkerEdgeColor','k','MarkerFaceColor','k');
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    legend_labels{end+1}=sprintf('%s, AUC=%0.2f ,norm pAUC=%0.3f', method, auc, norm_p_auc);
end
hold off;
legend(legend_labels);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
