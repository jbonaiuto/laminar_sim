function plot_roc_curves(subj_info, session_num, freq, snr, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10, ...
    'surf_dir', 'd:\pred_coding\surf');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
methods_to_plot=[1 4];
Nmeth=length(methodnames);

orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');

orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');

simmeshes={white_mesh,pial_mesh};
Nmesh=length(simmeshes);

labels=[-1*ones(params.nsims,1); 1*ones(params.nsims,1)];
f_scores=zeros(length(methodnames),params.nsims*Nmesh);
t_scores=zeros(length(methodnames),params.nsims*Nmesh);

figure();
hold on;
method_colors={'b','r'};
legend_labels={};
for i=1:length(methods_to_plot)
    methind=methods_to_plot(i);
    method=methodnames{methind};
    disp(method);
    
    data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),snr, params.dipole_moment));
    load(data_file);
    
    for simmeshind=1:Nmesh,    
        truotherF=squeeze(allcrossF(simmeshind,1:params.nsims,2,methind)-allcrossF(simmeshind,1:params.nsims,1,methind));
        f_scores(methind, (simmeshind-1)*params.nsims+1:simmeshind*params.nsims,:)=truotherF;
    end
    [x,y,t,auc]=compute_roc(f_scores(methind,:)',labels,params.nsims);
    plot(x,y,method_colors{i});
    if min(t)<0 && max(t)>0
        thresh_idx=find(t==0);
        h=plot(x(thresh_idx),y(thresh_idx),'o','MarkerEdgeColor','k','MarkerFaceColor','k');
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    legend_labels{end+1}=sprintf('F, %s, AUC=%0.2f',method, auc);        
end

data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id, num2str(session_num), sprintf('f%d_%d_SNR%d_dipolemoment%d', freq(1), freq(2), snr, params.dipole_moment));
for i=1:length(methods_to_plot),      
    methind=methods_to_plot(i);
    method=methodnames{methind};
    disp(method);
    
    labels=[-1*ones(params.nsims,1); 1*ones(params.nsims,1)];
    wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
        white_mesh, orig_pial_mesh, orig_white_mesh);
    t_scores(methind,:)=wmpial_t;
        
    [x,y,t,auc]=compute_roc(wmpial_t,labels,params.nsims);
    plot(x,y,sprintf('%s--',method_colors{i}));
    if min(t)<-0.01 && max(t)>0.01
        thresh_idx=find(t==0);
        h=plot(x(thresh_idx),y(thresh_idx),'o','MarkerEdgeColor','k','MarkerFaceColor','k');
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end 
    legend_labels{end+1}=sprintf('t, %s, AUC=%.2f', methodnames{methind}, auc);
end
hold off;
legend(legend_labels);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(sprintf('SNR=%ddB', snr));

disp('Comparing AUCs');

labels(find(labels<0))=0;
[pvalue Webb Wmsp]=pauc(f_scores(1,:)',f_scores(4,:)',labels);
disp(sprintf('Whole Brain, EBB-MSP, W_EBB=%.4f, W_MSP=%.4f, p=%.5f', Webb, Wmsp, pvalue));
[pvalue Webb Wmsp]=pauc(t_scores(1,:)',t_scores(4,:)',labels);
disp(sprintf('ROI, EBB-MSP, W_EBB=%.4f, W_MSP=%.4f, p=%.5f', Webb, Wmsp, pvalue));
[pvalue Wwb Wroi]=pauc(f_scores(1,:)',t_scores(1,:)',labels);
disp(sprintf('EBB, Whole brain-ROI, W_whole_brain=%.4f, W_roi=%.4f, p=%.5f', Wwb, Wroi, pvalue));
[pvalue Wwb Wroi]=pauc(f_scores(4,:)',t_scores(4,:)',labels);
disp(sprintf('MSP, Whole brain-ROI, W_whole_brain=%.4f, W_roi=%.4f, p=%.5f', Wwb, Wroi, pvalue));


    