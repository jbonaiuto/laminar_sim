function plot_ttest_paired_priors_results(subj_info, session_num, freq, snr, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10, ...
    'surf_dir', 'd:\pred_coding\surf', 'mri_dir', 'd:\pred_coding\mri');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
methods_to_plot=[4];
Nmeth=length(methodnames);

orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');
orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');
simmeshes={white_mesh,pial_mesh};
Nmesh=length(simmeshes);
data_dir=fullfile('D:\layer_sim\ttest_paired_prior_results', subj_info.subj_id, num2str(session_num), sprintf('f%d_%d_SNR%d_dipolemoment%d', freq(1), freq(2), snr, params.dipole_moment));

figure();
mesh_idx=[3 4 7 8;1 2 5 6];
%for methind=1:Nmeth,    
for i=1:length(methods_to_plot)
    methind=methods_to_plot(i);
    method=methodnames{methind};
    disp(method);
    
    wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
        white_mesh, orig_pial_mesh, orig_white_mesh);
    
    % For each simulated mesh
    for simmeshind=1:Nmesh,
        [path,file,ext]=fileparts(simmeshes{simmeshind});
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        simmeshname=y{2};
        
        subplot(4,2,mesh_idx(simmeshind,methind));
        hold on
        for simind=1:params.nsims
            color='b';
            if wmpial_t((simmeshind-1)*params.nsims+simind)<0
                color='r';
            end
            bar(simind,wmpial_t((simmeshind-1)*params.nsims+simind),color);
        end
        hold on
        plot([0 params.nsims+1], [1.69 1.69],'r--');
        plot([0 params.nsims+1], [-1.69 -1.69],'r--');
        xlim([0 params.nsims+1]);
        %ylim([-40 40]);
        xlabel('Simulation')
        ylabel('Mean t (pial-white)');
        title(sprintf('t Stat, %s, %s',methodnames{methind},simmeshname)); 
    end
end

figure();
hold all;
legend_labels={};
for i=1:length(methods_to_plot)
    methind=methods_to_plot(i);
    method=methodnames{methind};
    
    wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
        white_mesh, orig_pial_mesh, orig_white_mesh);
    
    labels=[-1*ones(params.nsims,1); 1*ones(params.nsims,1)];
    [x,y,t,auc]=perfcurve(labels,wmpial_t,'1');
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
    legend_labels{end+1}=sprintf('%s, AUC=%0.2f ,norm pAUC=%0.3f', method, auc, norm_p_auc);
end
hold off;
legend(legend_labels);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
