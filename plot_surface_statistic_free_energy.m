function plot_surface_statistic_free_energy(statistic, subj_info, session_num, freq, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'snr', 5, 'dipole_moment', 10,...
    'surf_dir', 'd:\pred_coding\surf', 'clip_vals', true, 'limits', [], ...
    'plot_dir', 'C:\Users\jbonai\Dropbox\meg\layer_sim');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

out_dir=fullfile(params.plot_dir,statistic);
mkdir(out_dir);

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
methods_to_plot=[1 4];
Nmeth=length(methodnames);

orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');

orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');

allmeshes=strvcat(white_mesh, pial_mesh);
Nmesh=size(allmeshes,1);

wm=gifti(white_mesh);
pial=gifti(pial_mesh);
wm_inflated=gifti(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed_inflated.surf.gii'));
pial_inflated=gifti(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed_inflated.surf.gii'));

pial_white_map=map_pial_to_white(white_mesh, pial_mesh, ...
        'mapType', 'link', 'origPial', orig_pial_mesh, ...
        'origWhite', orig_white_mesh);
    
switch statistic
    case 'thickness'
        [pial_statistic, wm_statistic]=compute_thickness(pial_mesh, ...
            white_mesh, orig_pial_mesh, orig_white_mesh);
    case 'curvature'
        pial_statistic=compute_curvature(pial_mesh, 'curvature_type', 'mean');
        wm_statistic=compute_curvature(white_mesh, 'curvature_type', 'mean');
    case 'depth'
        [pial_statistic,HS]=compute_sulcal_depth(pial);
        mapping=dsearchn(HS.vertices,wm.vertices);
        wm_statistic=sqrt(sum((wm.vertices-HS.vertices(mapping,:)).^2,2));
    case 'lead_field_norm'
        D=spm_eeg_load('C:\Users\jbonai\Dropbox\meg\layer_sim\rgb_1_pial.mat');
        pial_statistic=sqrt(sum(D.inv{1}.inverse.L.^2,1))';

        D=spm_eeg_load('C:\Users\jbonai\Dropbox\meg\layer_sim\rgb_1_white.mat');
        wm_statistic=sqrt(sum(D.inv{1}.inverse.L.^2,1))';
end

switch statistic
    case 'thickness'
        fig=figure();
        bins=linspace(min(pial_statistic),max(pial_statistic),100);
        [n,b] = histc(pial_statistic,bins);
        bar(bins,n./sum(n)*100.0,'histc');
        xlabel('Thickness');
        ylabel('% Vertices');
    otherwise
        fig=figure();
        plot(pial_statistic,wm_statistic(pial_white_map),'.');
        hold on
        min_val=min([pial_statistic(:) wm_statistic(:)]);
        max_val=max([pial_statistic(:) wm_statistic(:)]);
        plot([min_val max_val],[min_val max_val],'r-');
        mean((pial_statistic-wm_statistic(pial_white_map))./wm_statistic(pial_white_map))
        std((pial_statistic-wm_statistic(pial_white_map))./wm_statistic(pial_white_map))
        xlabel(sprintf('Pial %s',statistic));
        ylabel(sprintf('White matter %s',statistic));
end
%saveas(fig, fullfile(out_dir, sprintf('pial_white_%s.png', statistic)), 'png');
%figure2eps(fig, fullfile(out_dir, sprintf('pial_white_%s.eps', statistic)), 10, '-opengl');
        
[ax, metric_data]=plot_surface_metric(pial, pial_statistic, 'clip_vals', params.clip_vals, 'limits', params.limits);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);
fig=get(ax,'Parent');
%saveas(fig, fullfile(out_dir, sprintf('pial_%s.png', statistic)), 'png');

[ax, metric_data]=plot_surface_metric(pial_inflated, pial_statistic, 'clip_vals', params.clip_vals, 'limits', params.limits);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);
fig=get(ax,'Parent');
%saveas(fig, fullfile(out_dir, sprintf('pial_inflated_%s.png', statistic)), 'png');

[ax, metric_data]=plot_surface_metric(wm, wm_statistic, 'clip_vals', params.clip_vals, 'limits', params.limits);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);
fig=get(ax,'Parent');
%saveas(fig, fullfile(out_dir, sprintf('wm_%s.png', statistic)), 'png');

[ax, metric_data]=plot_surface_metric(wm_inflated, wm_statistic, 'clip_vals', params.clip_vals, 'limits', params.limits);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);
fig=get(ax,'Parent');
%saveas(fig, fullfile(out_dir, sprintf('wm_inflated_%s.png', statistic)), 'png');

nverts=size(wm.vertices,1);
rng(0);
simvertind=randperm(nverts); %% random list of vertex indices to simulate sources on

for i=1:length(methods_to_plot),       
    methind=methods_to_plot(i);
    method=methodnames{methind};
    fig=figure();
    hold on;
    
    allTrueOtherF=zeros(Nmesh*params.nsims,1);
    for simmeshind=1:Nmesh,    
        [path,file,ext]=fileparts(deblank(allmeshes(simmeshind,:)));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        
        data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),params.snr,params.dipole_moment));
        load(data_file);

        % F reconstructed on true - reconstructed on other
        % num simulations x number of folds
        truotherF=squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind)-allcrossF(simmeshind,1:params.nsims,(2-simmeshind)+1,methind));
        allTrueOtherF((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=truotherF;
        
        switch simmeshind
            case 1
                sim_stats=wm_statistic(simvertind(1:params.nsims));
                color='r';
            case 2
                sim_stats=pial_statistic(simvertind(1:params.nsims));
                color='b';
        end
        plot(sim_stats,truotherF,'o','MarkerEdgeColor','none','MarkerFaceColor',color);        
    end
    
    all_stat=[wm_statistic(simvertind(1:params.nsims)); pial_statistic(simvertind(1:params.nsims))];
    [rho,pval] = corr(all_stat, allTrueOtherF,'type','Spearman');
    pPoly = polyfit(all_stat, allTrueOtherF, 1); % Linear fit of xdata vs ydata
    linePointsX = [min(all_stat) max(all_stat)]; % find left and right x values
    linePointsY = polyval(pPoly,linePointsX); 
    plot (linePointsX,linePointsY,'--k')
    
    hold off;
    xlabel(statistic)
    ylabel('Free energy diff (correct-incorrect)');
    title(sprintf('Free energy, %s, Rho=%.5f, p=%.5f',methodnames{methind}, rho, pval));        
    
    %saveas(fig, fullfile(out_dir, sprintf('%s_free_energy_%s_%d-%dHz_%ddb.png', statistic, method, freq(1), freq(2), params.snr)), 'png');
    %figure2eps(fig, fullfile(out_dir, sprintf('%s_free_energy_%s_%d-%dHz_%ddb.eps', statistic, method, freq(1), freq(2), params.snr)), 10, '-opengl');
end