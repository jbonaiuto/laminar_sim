function tstat=compute_roi_t_stat(file_prefix, pial_meshname, wm_meshname, varargin)

% Parse inputs
defaults = struct('mapType', 'link', 'recompute', false, 'origPial', '',...
    'origWhite', '', 'threshold',[]);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

[filepath filename ext]=fileparts(file_prefix);
% Files containing t-statistics and pial-wm diff
pial_t_filename=fullfile(filepath, sprintf('pial.%s.t.gii', filename));
wm_t_filename=fullfile(filepath, sprintf('white.%s.t.gii', filename));
pial_wm_diff_filename=fullfile(filepath, sprintf('pial-white.%s.gii', filename));
pial_white_map=map_pial_to_white(wm_meshname, pial_meshname, ...
        'mapType', params.mapType, 'origPial', params.origPial, ...
        'origWhite', params.origWhite);
    
% If tstat files or pial_wm file do not exist or recomputing
if exist(pial_t_filename,'file')~=2 || exist(wm_t_filename,'file')~=2 || exist(pial_wm_diff_filename, 'file')~=2 || params.recompute
    fprintf('Recomputing t-stat for %s\n', file_prefix);
    
    % Load pial and white matter surface (woi-baseline)
    pial_diff=gifti(fullfile(filepath,sprintf('pial.%s.gii',filename)));
    wm_diff=gifti(fullfile(filepath,sprintf('white.%s.gii',filename)));
            
    % Compute correction for pial surface t-test
    correction = 100.0*max(nanstd(pial_diff.cdata(:,:),[],2)); 
    % Run pial surface t-test
    [tstat,pvals]=ttest_corrected(pial_diff.cdata(:,:)','correction',correction);    
    pial_tvals=tstat';
    write_metric_gifti(pial_t_filename, pial_tvals);
    
    % Compute correction for wm surface t-test
    correction = 100.0*max(nanstd(wm_diff.cdata(:,:),[],2));
    % Run wm surface t-test
    [tstat,pvals]=ttest_corrected(wm_diff.cdata(:,:)','correction',correction);
    wm_tvals=tstat';
    write_metric_gifti(wm_t_filename, wm_tvals);
            
    % Compute pial-white difference
    pial_wm_diff=abs(pial_diff.cdata(:,:))-abs(wm_diff.cdata(pial_white_map,:));
    write_metric_gifti(pial_wm_diff_filename, pial_wm_diff);
else % Otherwise load data from files
    x=gifti(pial_t_filename);
    pial_tvals=x.cdata(:);
    x=gifti(wm_t_filename);
    wm_tvals=x.cdata(:);
    x=gifti(pial_wm_diff_filename);
    pial_wm_diff=x.cdata(:,:);
end

% If threshold not specified - use 10% of max-min
if length(params.threshold)==0
    min_t=min([pial_tvals; wm_tvals]);
    max_t=max([pial_tvals; wm_tvals]);
    params.threshold=min_t+0.1*(max_t-min_t);
end
            
% Create mask
[pial_mask,wm_mask,mask]=get_pial_wm_mask(pial_tvals, wm_tvals, ...
    params.threshold, pial_white_map);
fprintf('ROI size=%d\n',length(mask));
% Get mean pial-wm in ROI
pial_wm_roi_diff=mean(pial_wm_diff(mask,:));
% Compute t-stat correction
correction=100*max(nanstd(pial_wm_diff,[],2));
% Perform ROI t-stat
[tstat,pvals]=ttest_corrected(pial_wm_roi_diff','correction',correction);
