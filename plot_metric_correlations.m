function plot_metric_correlations(subj_info, session_num, freq, snr,...
    method_idx, varargin)
% PLOT_METRIC_CORRELATIONS  Plot correlations between t-statistics, free
% energy, and CV error
%
% Use as
%   plot_metric_correlations(subj_info, 1, [10 30], -20, 1)
% where the first argument is the subject info structure (from create_subjects),
% the second is the session numner, the third is the frequency range of the
% simulated data (Hz), the fourth is the SNR (db), and the fifth is the
% method index (1=EBB, 2=IID, 3=COH, 4=MSP)
% 
%   plot_metric_correlations(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * sources - 'both' (default) or ('pial','white') - which sources to
%    plot
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces

% Parse inputs
defaults = struct('sources','both','nsims', 60, 'dipole_moment', 10,...
    'surf_dir', 'd:\pred_coding\surf');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames={'EBB','IID','COH','MSP'};
method=methodnames{method_idx};
cv_method_idx=method_idx;
% Only ran EBB and MSP for cross validation error
if method_idx==4
    cv_method_idx=2;
end

% Original and downsampled white matter surface
orig_white_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_white.hires.deformed.surf.gii');

% Original and downsampled pial surface
orig_pial_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_pial.hires.deformed.surf.gii');

% Load cross-validation results
fname=sprintf('allcrossErr_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),...
    snr,params.dipole_moment);
crosserr_file=fullfile('d:/layer_sim/results/',subj_info.subj_id,num2str(session_num),...
    fname);
load(crosserr_file);
% Normalize and average of left-out channels
allCrossErr=squeeze(mean(results.allCrossErr./results.allRMS,6));
% Average over folds
meancrossErr=squeeze(mean(allCrossErr,5));

% Load free energy results
fname=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),...
    snr,params.dipole_moment);
f_file=fullfile('d:/layer_sim/results/',subj_info.subj_id,num2str(session_num),...
    fname);
load(f_file);

% Load ROI t-statistics
fname=sprintf('f%d_%d_SNR%d_dipolemoment%d',freq(1),freq(2),snr,...
    params.dipole_moment);
data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id, num2str(session_num), ...
    fname);
wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
        white_mesh, orig_pial_mesh, orig_white_mesh);
        
% Get white matter and pial simulation t-stats
wmT=wmpial_t(1:params.nsims);
pialT=wmpial_t(params.nsims+1:2*params.nsims);
% Get CV error difference: White - pial (because small CV err is better, large F value is better
wmmeancrossErr=squeeze(meancrossErr(1,1:params.nsims,1,cv_method_idx)-meancrossErr(1,1:params.nsims,2,cv_method_idx)).*100;
pialmeancrossErr=squeeze(meancrossErr(2,1:params.nsims,1,cv_method_idx)-meancrossErr(2,1:params.nsims,2,cv_method_idx)).*100;
% Get F difference: Pial-white
wmF=squeeze(allcrossF(1,1:params.nsims,2,method_idx)-allcrossF(1,1:params.nsims,1,method_idx));
pialF=squeeze(allcrossF(2,1:params.nsims,2,method_idx)-allcrossF(2,1:params.nsims,1,method_idx));

% Show both surface simulations, just pial or just white
switch params.sources
    case 'both'
        f=[wmF pialF]';
        cverr=[wmmeancrossErr pialmeancrossErr]';
    case 'white'
        f=wmF';
        cverr=wmmeancrossErr';
        wmpial_t=wmT;
    case 'pial'
        f=pialF';
        cverr=pialmeancrossErr';
        wmpial_t=pialT;
end

figure();
% Free energy, CV err, T
subplot(2,2,1); % CV err x free energy
hold on;
if strcmp(params.sources,'both') || strcmp(params.sources,'white')
    plot(wmmeancrossErr,wmF,'o','MarkerEdgeColor','none','MarkerFaceColor','r');
end
if strcmp(params.sources,'both') || strcmp(params.sources,'pial')
    plot(pialmeancrossErr,pialF,'o','MarkerEdgeColor','none','MarkerFaceColor','b');
end
[rho,pval] = corr(cverr, f,'type','Spearman');    
pPoly = polyfit(cverr, f, 1); % Linear fit of xdata vs ydata
linePointsX = [min(cverr) max(cverr)]; % find left and right x values
linePointsY = polyval(pPoly,[min(cverr),max(cverr)]); 
plot (linePointsX,linePointsY,'--k')
hold off;
xlabel('CV Err');
ylabel('F');
title(sprintf('%s, rho=%.2f, p=%0.5f\n', method,rho,pval));
legend({'Deep source','Superficial source'});

subplot(2,2,2); % t x free energy
hold on;
switch method_idx
    % Split MSP because there are two clusters
    case 4
        plot(wmT,wmF,'o','MarkerEdgeColor','none','MarkerFaceColor','r');
        [rhow,pvalw] = corr(wmT, wmF','type','Spearman');    
        pPoly = polyfit(wmT, wmF', 1); % Linear fit of xdata vs ydata
        linePointsX = [min(wmT) max(wmT)]; % find left and right x values
        linePointsY = polyval(pPoly,[min(wmT),max(wmT)]); 
        plot (linePointsX,linePointsY,'--k')
        
        plot(pialT,pialF,'o','MarkerEdgeColor','none','MarkerFaceColor','b');
        [rhop,pvalp] = corr(pialT, pialF','type','Spearman');    
        pPoly = polyfit(pialT, pialF', 1); % Linear fit of xdata vs ydata
        linePointsX = [min(pialT) max(pialT)]; % find left and right x values
        linePointsY = polyval(pPoly,[min(pialT),max(pialT)]); 
        plot (linePointsX,linePointsY,'--k')
        title(sprintf('%s, pial: rho=%.2f, p=%0.5f, white: rho=%.2f, p=%0.5f\n', method,rhop,pvalp,rhow,pvalw));
    otherwise
        if strcmp(params.sources,'both') || strcmp(params.sources,'white')
            plot(wmT,wmF,'o','MarkerEdgeColor','none','MarkerFaceColor','r');
        end
        if strcmp(params.sources,'both') || strcmp(params.sources,'pial')
            plot(pialT,pialF,'o','MarkerEdgeColor','none','MarkerFaceColor','b');
        end
        [rho,pval] = corr(wmpial_t, f,'type','Spearman');    
        pPoly = polyfit(wmpial_t, f, 1); % Linear fit of xdata vs ydata
        linePointsX = [min(wmpial_t) max(wmpial_t)]; % find left and right x values
        linePointsY = polyval(pPoly,[min(wmpial_t),max(wmpial_t)]); 
        plot (linePointsX,linePointsY,'--k')
        title(sprintf('%s, rho=%.2f, p=%0.5f\n', method,rho,pval));
end
hold off;
xlabel('t');
ylabel('F');


subplot(2,2,4); % CV x t
hold on;
switch method_idx
    % Split MSP because there are two clusters
    case 4
        plot(wmT,wmmeancrossErr,'o','MarkerEdgeColor','none','MarkerFaceColor','r');
        [rhow,pvalw] = corr(wmT, wmmeancrossErr','type','Spearman');    
        pPoly = polyfit(wmT, wmmeancrossErr', 1); % Linear fit of xdata vs ydata
        linePointsX = [min(wmT) max(wmT)]; % find left and right x values
        linePointsY = polyval(pPoly,[min(wmT),max(wmT)]); 
        plot (linePointsX,linePointsY,'--k');
        
        plot(pialT,pialmeancrossErr,'o','MarkerEdgeColor','none','MarkerFaceColor','b');
        [rhop,pvalp] = corr(pialT, pialmeancrossErr','type','Spearman');    
        pPoly = polyfit(pialT, pialmeancrossErr', 1); % Linear fit of xdata vs ydata
        linePointsX = [min(pialT) max(pialT)]; % find left and right x values
        linePointsY = polyval(pPoly,[min(pialT),max(pialT)]); 
        plot (linePointsX,linePointsY,'--k');
        title(sprintf('%s, pial: rho=%.2f, p=%0.5f, white: rho=%.2f, p=%0.5f\n', method,rhop,pvalp,rhow,pvalw));
    otherwise
        if strcmp(params.sources,'both') || strcmp(params.sources,'white')
            plot(wmT,wmmeancrossErr,'o','MarkerEdgeColor','none','MarkerFaceColor','r');
        end
        if strcmp(params.sources,'both') || strcmp(params.sources,'pial')
            plot(pialT,pialmeancrossErr,'o','MarkerEdgeColor','none','MarkerFaceColor','b');
        end
        [rho,pval] = corr(wmpial_t, cverr,'type','Spearman');    
        pPoly = polyfit(wmpial_t, cverr, 1); % Linear fit of xdata vs ydata
        linePointsX = [min(wmpial_t) max(wmpial_t)]; % find left and right x values
        linePointsY = polyval(pPoly,[min(wmpial_t),max(wmpial_t)]); 
        plot (linePointsX,linePointsY,'--k');
        title(sprintf('%s, rho=%.2f, p=%0.5f\n', method,rho,pval));
end
hold off;
xlabel('t');
ylabel('CV Err');

