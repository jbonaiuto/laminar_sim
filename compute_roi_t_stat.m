function tstat=compute_roi_t_stat(file_prefix, pial_meshname, wm_meshname, varargin)

% Parse inputs
defaults = struct('mapType', 'link', 'recompute', false, 'origPial', '',...
    'origWhite', '');  %define default values
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
            
    % Run pial surface t-test
    varpop=nanvar([pial_diff.cdata(:,:) wm_diff.cdata(:,:)],[],2);
    %varpop=nanvar(pial_diff.cdata(:,:),[],2);
    [tstat,pvals]=ttest_corrected(pial_diff.cdata(:,:)','correction',.01*max(varpop));
    %[tstat,pvals]=ttest_corrected(pial_diff.cdata(:,:)');
    pial_tvals=tstat';
    %[H,pvals,ci,STATS]=ttest(pial_diff.cdata(:,:)');
    %pial_tvals=STATS.tstat';
    write_metric_gifti(pial_t_filename, pial_tvals);
    
    % Run wm surface t-test
    %varpop=nanvar(wm_diff.cdata(:,:),[],2);
    [tstat,pvals]=ttest_corrected(wm_diff.cdata(:,:)','correction',.01*max(varpop));
    %[tstat,pvals]=ttest_corrected(wm_diff.cdata(:,:)');
    wm_tvals=tstat';
    %[H,pvals,ci,STATS]=ttest(wm_diff.cdata(:,:)');
    %wm_tvals=STATS.tstat';
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

pial_threshold=prctile(pial_tvals(~isinf(pial_tvals)),95);
% Create pial and white masks and mapped white mask
pial_mask=find(pial_tvals>pial_threshold & ~isinf(pial_tvals));

wm_threshold=prctile(wm_tvals(~isinf(wm_tvals)),95);
wm_mask=find(wm_tvals>wm_threshold & ~isinf(wm_tvals));
mapped_wm_tvals=wm_tvals(pial_white_map);
mapped_wm_mask=find(mapped_wm_tvals>wm_threshold & ~isinf(mapped_wm_tvals));
mask=union(pial_mask, mapped_wm_mask);
        
% Get mean pial-wm in ROI
pial_wm_roi_diff=mean(pial_wm_diff(mask,:));
% Perform ROI t-stat
%varpop=nanvar([pial_diff.cdata(:,:) wm_diff.cdata(:,:)],[],2);
%[tstat,pvals]=ttest_corrected(pial_wm_roi_diff','correction',10*max(varpop));
[H,pvals,ci,STATS]=ttest(pial_wm_roi_diff');
tstat=STATS.tstat;
fprintf('ROI size=%d, tstat=%.2f\n',length(mask),tstat);
