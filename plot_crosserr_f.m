function plot_crosserr_f(subj_info, session_num, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'surf_dir', 'd:\pred_coding\surf', 'mri_dir', 'd:\pred_coding\mri');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

crosserr_file=fullfile('d:/layer_sim/results/',subj_info.subj_id,num2str(session_num),'allcrossErr_f15_30_SNR5.mat');
load(crosserr_file);
f_file=fullfile('d:/layer_sim/results/',subj_info.subj_id,num2str(session_num),'allcrossF_f15_30_SNR5_dipolemoment10.mat');
load(f_file);

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now

white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');
simmeshes={white_mesh,pial_mesh};


meancrossErr=squeeze(mean(allcrossErr,5));
% White - pial (because small CV err is better, large F value is better
wmmeancrossErr=squeeze(meancrossErr(1,:,1,:)-meancrossErr(1,:,2,:));
pialmeancrossErr=squeeze(meancrossErr(2,:,1,:)-meancrossErr(2,:,2,:));
wmF=squeeze(allcrossF(1,:,2,:)-allcrossF(1,:,1,:));
pialF=squeeze(allcrossF(2,:,2,:)-allcrossF(2,:,1,:));

for methind=1:length(methodnames)
    figure();
    cve=[wmmeancrossErr(:,methind); pialmeancrossErr(:,methind)];
    f=[wmF(:,methind); pialF(:,methind)];

    hold on;
    plot(wmmeancrossErr(:,methind),wmF(:,methind),'or');
    plot(pialmeancrossErr(:,methind),pialF(:,methind),'ob');

    [rho,pval] = corr(cve, f,'type','Spearman');    
    pPoly = polyfit(cve, f, 1); % Linear fit of xdata vs ydata
    linePointsX = [min(cve) max(cve)]; % find left and right x values
    linePointsY = polyval(pPoly,[min(cve),max(cve)]); 
    plot (linePointsX,linePointsY,'--k')
    hold off;
    xlabel('Cross Val Err');
    ylabel('F');
    title(sprintf('%s, rho=%.2f, p=%0.5f\n', methodnames{methind},rho,pval));
    legend({'Deep source','Superficial source'});
end

