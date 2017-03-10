function plot_perc_correct(subj_info, session_num, freq, snr, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10, 'surf_dir', 'd:\pred_coding\surf', ...
    'mri_dir', 'd:\pred_coding\mri','threshold',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
Nmeth=length(methodnames);

orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');

orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');

simmeshes={white_mesh,pial_mesh};
Nmesh=length(simmeshes);

dof=514;
alpha=1.0-0.05/2;
t_thresh=tinv(alpha, dof);

perc_correct=zeros(2,Nmeth);
perc_correct_stderr=zeros(2,Nmeth);
for methind=1:Nmeth,           
    data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),snr, params.dipole_moment));
    load(data_file);
    
    correct=zeros(1,Nmesh*params.nsims);
    for simmeshind=1:Nmesh,    
        truotherF=squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind)-allcrossF(simmeshind,1:params.nsims,(2-simmeshind)+1,methind));
        if params.threshold
            correct((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=truotherF>3;
        else
            correct((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=truotherF>0;
        end
    end
    perc_correct(1,methind)=mean(correct);
    perc_correct_stderr(1,methind)=std(correct)/sqrt(length(correct));
end

data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id, num2str(session_num), sprintf('f%d_%d_SNR%d_dipolemoment%d', freq(1), freq(2), snr, params.dipole_moment));
for methind=1:Nmeth,      
    method=methodnames{methind};
    correct=zeros(1,Nmesh*params.nsims);
    
    wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
        white_mesh, orig_pial_mesh, orig_white_mesh);
    
    for simmeshind=1:Nmesh,    
        for s=1:params.nsims
            tstat=wmpial_t((simmeshind-1)*params.nsims+s);
            if params.threshold
                if (simmeshind==1 && tstat<-t_thresh) || (simmeshind==2 && tstat>t_thresh)
                    correct((simmeshind-1)*params.nsims+s)=1;
                end
            else
                if (simmeshind==1 && tstat<0) || (simmeshind==2 && tstat>0)
                    correct((simmeshind-1)*params.nsims+s)=1;
                end
            end
        end
    end    
    perc_correct(2,methind)=mean(correct);
    perc_correct_stderr(2,methind)=std(correct)/sqrt(length(correct));
end

figure();
barwitherr(perc_correct_stderr,perc_correct);
set(gca,'XTickLabel',{'Free Energy','t Test'});
legend(methodnames);
ylim([0 1]);
ylabel('% Correct');
