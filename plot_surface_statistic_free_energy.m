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

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
methods_to_plot=[1 4];
Nmeth=length(methodnames);

[pial_statistic,wm_statistic]=plot_surface_statistic(statistic, subj_info,...
    'surf_dir', params.surf_dir, 'clip_vals', params.clip_vals, ...
    'limits', params.limits, 'plot_dir', params.plot_dir);

nverts=size(pial_statistic,1);
rng(0);
simvertind=randperm(nverts); %% random list of vertex indices to simulate sources on
Nmesh=2;

for i=1:length(methods_to_plot),       
    methind=methods_to_plot(i);
    method=methodnames{methind};
    fig=figure();
    hold on;
    
    allTrueOtherF=zeros(Nmesh*params.nsims,1);
    for simmeshind=1:Nmesh,    
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
    
    saveas(fig, fullfile(params.plot_dir, sprintf('%s_free_energy_%s_%d-%dHz_%ddb.png', statistic, method, freq(1), freq(2), params.snr)), 'png');
    figure2eps(fig, fullfile(params.plot_dir, sprintf('%s_free_energy_%s_%d-%dHz_%ddb.eps', statistic, method, freq(1), freq(2), params.snr)), 10, '-opengl');
end
