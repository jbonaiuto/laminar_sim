function plot_surface_statistic_free_energy(statistic, subj_info,...
    session_num, freq, varargin)
% PLOT_SURFACE_STATISTIC_FREE_ENERGY  Plot relationship between surface statistic
% and free energy
%
% Use as
%   plot_surface_statistic_free_energy('depth', subjects(1), 1, [10 30])
% where the first argument is the statistic (thickness, depth, curvature, or lead_field_norm,
% the second is the subject info structure (from create_subjects), the third is the
% session number, and the fourth is the frequency range.
% 
%   plot_surface_statistic_free_energy(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * snr - -20 (default) or integer - SNR (db)
%    * dipole_moment - 10 (default) or interger - moment of simulated
%    dipole
%    * surf_dir - directory containing subject surfaces
%    * clip_vals - true (default) or boolean - whether or not to clip
%    values at 95% positive and negative
%    * limits - [] (default) or vector - color map limits
%    * plot_dir - directory to save plots

% Parse inputs
defaults = struct('nsims', 60, 'snr', -20, 'dipole_moment', 10,...
    'surf_dir', 'd:\pred_coding\surf', 'clip_vals', true, 'limits', [], ...
    'plot_dir', 'C:\Users\jbonai\Dropbox\meg\layer_sim');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames={'EBB','IID','COH','MSP'};

% Get surface statistic for each surface
[pial_statistic,wm_statistic]=plot_surface_statistic(statistic, subj_info,...
    'surf_dir', params.surf_dir, 'clip_vals', params.clip_vals, ...
    'limits', params.limits, 'plot_dir', params.plot_dir);

nverts=size(pial_statistic,1);
rng(0);
simvertind=randperm(nverts); %% random list of vertex indices to simulate sources on
Nmesh=2;

% Plot EBB and MSP
methods_to_plot=[1 4];
for i=1:length(methods_to_plot),       
    methind=methods_to_plot(i);
    figure();
    hold on;
    
    allTrueOtherF=zeros(Nmesh*params.nsims,1);
    % Load free energy results
    fname=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),...
        freq(2),params.snr,params.dipole_moment);
    data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, ...
        num2str(session_num), fname);
    load(data_file);

    for simmeshind=1:Nmesh,            
        % F reconstructed on true - reconstructed on other
        trueF=squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind));
        otherF=squeeze(allcrossF(simmeshind,1:params.nsims,(2-simmeshind)+1,methind));
        allTrueOtherF((simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=trueF-otherF;
        
        switch simmeshind
            case 1
                sim_stats=wm_statistic(simvertind(1:params.nsims));
                color='r';
            case 2
                sim_stats=pial_statistic(simvertind(1:params.nsims));
                color='b';
        end
        plot(sim_stats,trueF-otherF,'o','MarkerEdgeColor','none','MarkerFaceColor',color);        
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
end
