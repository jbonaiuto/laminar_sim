function plot_snr_perc_correct(subj_info, session_num, freq, snrs, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10, 'surf_dir', 'd:\pred_coding\surf', ...
    'mri_dir', 'd:\pred_coding\mri','threshold',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

dof=514;
alpha=1.0-0.05/2;
t_thresh=tinv(alpha, dof);

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
methods_to_plot=[1 4];
Nmeth=length(methodnames);

orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');

orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');

simmeshes={white_mesh,pial_mesh};
Nmesh=length(simmeshes);

figure();
hold on;
%method_colors={'b','r','g','c'};
method_colors={'b','r'};
legend_labels={};
%for methind=1:length(methodnames)
for i=1:length(methods_to_plot)
    methind=methods_to_plot(i);
    method=methodnames{methind};
    
    perc_correct=[];
    stderr_perc_correct=[];
    snrs_to_use=[];
    
    for s=1:length(snrs)
        snr=snrs(s);
        data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),snr, params.dipole_moment));
        load(data_file);

        correct=[];
        for simmeshind=1:Nmesh,    
            truotherF=squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind)-allcrossF(simmeshind,1:params.nsims,(2-simmeshind)+1,methind));
            if params.threshold
                for x=1:params.nsims
                    if abs(truotherF(x))>3
                        correct(end+1)=truotherF(x)>3;
                    end
                end
            else
                correct(end+1:end+length(truotherF))=truotherF>0;
            end
        end
        if length(correct)>.05*(Nmesh*params.nsims)
            perc_correct(end+1)=mean(correct);
            stderr_perc_correct(end+1)=std(correct)/sqrt(length(correct));
            snrs_to_use(end+1)=snr;
        end
    end
    errorbar(snrs_to_use,perc_correct.*100,stderr_perc_correct.*100,method_colors{i});
    legend_labels{end+1}=sprintf('F, %s',method);        
end

for i=1:length(methods_to_plot)
    methind=methods_to_plot(i);
    method=methodnames{methind};
    
    perc_correct=[];
    stderr_perc_correct=[];
    snrs_to_use=[];
    
    for j=1:length(snrs)
        snr=snrs(j);
        data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id, num2str(session_num), sprintf('f%d_%d_SNR%d_dipolemoment%d', freq(1), freq(2), snr, params.dipole_moment));
        correct=zeros(1,Nmesh*params.nsims);
        
        wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
            white_mesh, orig_pial_mesh, orig_white_mesh);
        
        correct=[];
        
        for simmeshind=1:Nmesh,    
            for s=1:params.nsims
                tstat=wmpial_t((simmeshind-1)*params.nsims+s);
                if params.threshold
                    if abs(tstat)>t_thresh
                        correct(end+1)=(simmeshind==1 && tstat<-t_thresh) || (simmeshind==2 && tstat>t_thresh);
                    end
                else
                    correct(end+1)=(simmeshind==1 && tstat<0) || (simmeshind==2 && tstat>0);
                end
            end
        end  
        if length(correct)>.05*(Nmesh*params.nsims)
            perc_correct(end+1)=mean(correct);
            stderr_perc_correct(end+1)=std(correct)/sqrt(length(correct));
            snrs_to_use(end+1)=snr;
        end
    end
    errorbar(snrs_to_use,perc_correct.*100,stderr_perc_correct.*100,sprintf('%s--',method_colors{i}));
    legend_labels{end+1}=sprintf('t, %s',method);
end

hold off;
legend(legend_labels);
xlabel('SNR');
ylabel('% Correct');
%ylim([50 100]);
ylim([40 105]);