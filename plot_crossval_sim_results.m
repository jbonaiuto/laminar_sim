function plot_crossval_sim_results(subj_info, session_num, freq, snrs, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'surf_dir', 'd:\pred_coding\surf', 'mri_dir', 'd:\pred_coding\mri');  %define default values
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

figure();
mesh_idx=[3 4 7 8;1 2 5 6];
for methind=1:Nmeth,       
    for simmeshind=1:Nmesh,    
        subplot(4,2,mesh_idx(simmeshind,methind));
        [path,file,ext]=fileparts(deblank(allmeshes(simmeshind,:)));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        simmeshname=y{2};
        
        allTruOtherErr=[];
        for s=1:length(snrs)
            snr=snrs(s);
            data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossErr_f%d_%d_SNR%d.mat',freq(1),freq(2),snr));
            load(data_file);

            % F reconstructed on true - reconstructed on other
            % num simulations x number of folds
            truotherErr=squeeze(mean(allcrossErr(simmeshind,1:params.nsims,1,methind,:),5)-mean(allcrossErr(simmeshind,1:params.nsims,2,methind,:),5));
            allTruOtherErr(end+1,:)=truotherErr;                        
        end
        bar(allTruOtherErr')
        hold on
        %plot([0 params.nsims+1],[3 3],'r--');
        %plot([0 params.nsims+1],[-3 -3],'r--');
        xlim([0 params.nsims+1]);
        if simmeshind==2
            legend('SNR=-5dB','SNR=0dB','SNR=5dB');
        end
        xlabel('Simulation')
        ylabel('Free energy diff (pial-white)');
        title(sprintf('Free energy, %s, %s',methodnames{methind},simmeshname));        
    end

end

for s=1:length(snrs)
    snr=snrs(s);
    
    figure();
    hold all;
    legend_labels={};
    
    for methind=1:Nmeth,           
        method=methodnames{methind};
        
        data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossErr_f%d_%d_SNR%d.mat',freq(1),freq(2),snr));
        load(data_file);
        labels=[-1*ones(params.nsims,1); 1*ones(params.nsims,1)];
        scores=zeros(params.nsims*Nmesh,1);    
        for simmeshind=1:Nmesh,    
            truotherErr=squeeze(mean(allcrossErr(simmeshind,1:params.nsims,2,methind,:),5)-mean(allcrossErr(simmeshind,1:params.nsims,1,methind,:)));
            scores((simmeshind-1)*params.nsims+1:simmeshind*params.nsims,:)=truotherErr;
        end
        [x,y,t,auc]=perfcurve(labels,scores,'1');
        plot(x,y);
        legend_labels{end+1}=sprintf('%s, AUC=%0.2f',method, auc);        
    end
    hold off;
    legend(legend_labels);
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    title(sprintf('SNR=%ddB', snr));
end
