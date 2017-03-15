function plot_perc_classified(subj_info, session_num, freq, snr, varargin)

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

perc_pial=zeros(2,Nmeth);
perc_pial_stderr=zeros(2,Nmeth);
perc_white=zeros(2,Nmeth);
perc_white_stderr=zeros(2,Nmeth);
for methind=1:Nmeth,           
    data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),snr, params.dipole_moment));
    load(data_file);
    
    pial=zeros(1,Nmesh*params.nsims);
    white=zeros(1,Nmesh*params.nsims);
    pialwhiteF=squeeze(allcrossF(:,1:params.nsims,2,methind)-allcrossF(:,1:params.nsims,1,methind));
    if params.threshold
        pial=pialwhiteF>3;
        white=pialwhiteF<-3;
    else
        pial=pialwhiteF>0;
        white=pialwhiteF<0;
    end
    perc_pial(1,methind)=mean(pial(:));
    perc_pial_stderr(1,methind)=std(pial(:))/sqrt(length(pial(:)));
    perc_white(1,methind)=mean(white(:));
    perc_white_stderr(1,methind)=std(white(:))/sqrt(length(white(:)));
end

data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id, num2str(session_num), sprintf('f%d_%d_SNR%d_dipolemoment%d', freq(1), freq(2), snr, params.dipole_moment));
for methind=1:Nmeth,      
    method=methodnames{methind};
    pial=zeros(1,Nmesh*params.nsims);
    white=zeros(1,Nmesh*params.nsims);
    
    wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
        white_mesh, orig_pial_mesh, orig_white_mesh);
    
    for simmeshind=1:Nmesh,    
        for s=1:params.nsims
            tstat=wmpial_t((simmeshind-1)*params.nsims+s);
            if params.threshold
                if tstat>t_thresh
                    pial((simmeshind-1)*params.nsims+s)=1;
                end
                if tstat<-t_thresh
                    white((simmeshind-1)*params.nsims+s)=1;
                end
            else
                if tstat>0
                    pial((simmeshind-1)*params.nsims+s)=1;
                end
                if tstat<0
                    white((simmeshind-1)*params.nsims+s)=1;
                end
            end
        end
    end             
    perc_pial(2,methind)=mean(pial);
    perc_pial_stderr(2,methind)=std(pial)/sqrt(length(pial));
    perc_white(2,methind)=mean(white);
    perc_white_stderr(2,methind)=std(white)/sqrt(length(white));
end

figure();
subplot(2,1,1);
barwitherr(perc_pial_stderr.*100,perc_pial.*100);
set(gca,'XTickLabel',{'Free Energy','t Test'});
legend(methodnames);
ylim([0 100]);
ylabel('% Pial');

subplot(2,1,2);
barwitherr(perc_white_stderr.*100,perc_white.*100);
set(gca,'XTickLabel',{'Free Energy','t Test'});
legend(methodnames);
ylim([0 100]);
ylabel('% White');
