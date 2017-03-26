function plot_crosserr_tval(subj_info, session_num, freq, snr, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10, 'surf_dir', 'd:\pred_coding\surf');  %define default values
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

crosserr_file=fullfile('d:/layer_sim/results/',subj_info.subj_id,num2str(session_num),...
    sprintf('allcrossErr_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),snr,params.dipole_moment));
load(crosserr_file);
    
meancrossErr=squeeze(mean(allcrossErr,5));
% White - pial (because small CV err is better, large F value is better
wmmeancrossErr=squeeze(meancrossErr(1,:,1,:)-meancrossErr(1,:,2,:));
pialmeancrossErr=squeeze(meancrossErr(2,:,1,:)-meancrossErr(2,:,2,:));

data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id, num2str(session_num), ...
    sprintf('f%d_%d_SNR%d_dipolemoment%d',freq(1),freq(2),snr,params.dipole_moment));

for methind=1:Nmeth,           
    method=methodnames{methind};
    figure();    
    
    wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
        white_mesh, orig_pial_mesh, orig_white_mesh);
           
    hold on;
    plot(wmmeancrossErr(:,methind),wmpial_t(1:params.nsims),'or');
    plot(pialmeancrossErr(:,methind),wmpial_t(params.nsims+1:end),'ob');

    cve=[wmmeancrossErr(:,methind); pialmeancrossErr(:,methind)];
    [rho,pval] = corr(cve, wmpial_t,'type','Spearman');    
    pPoly = polyfit(cve, wmpial_t, 1); % Linear fit of xdata vs ydata
    linePointsX = [min(cve) max(cve)]; % find left and right x values
    linePointsY = polyval(pPoly,[min(cve),max(cve)]); 
    plot (linePointsX,linePointsY,'--k')
    hold off;
    xlabel('Cross validation error');
    ylabel('t');
    title(sprintf('%s, rho=%.2f, p=%0.5f\n', methodnames{methind},rho,pval));
    legend({'Deep source','Superficial source'});
end
