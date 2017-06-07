function wmpial_t=get_wmpial_t(data_dir, method, nsims, pial_mesh, ...
    white_mesh, orig_pial_mesh, orig_white_mesh, varargin)

% Parse inputs
defaults = struct('recompute',false, 'recompute_trials', false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

file_name=fullfile(data_dir, sprintf('trials.%s.sim_mesh.mat', method));
if exist(file_name,'file')~=2 || params.recompute || params.recompute_trials
    wmpial_t=zeros(2*nsims,1);
    for simmeshind=1:2,    
        for s=1:nsims
            file_prefix=fullfile(data_dir,sprintf('%s.sim_mesh%d_source%d.gii',method,simmeshind,s));
            tstat=compute_roi_t_stat(file_prefix, pial_mesh, white_mesh,...
                'mapType', 'link', 'origPial', orig_pial_mesh,...
                'origWhite', orig_white_mesh,...
                'recompute', params.recompute_trials); 
            wmpial_t((simmeshind-1)*nsims+s)=tstat;
        end
    end
    save(file_name, 'wmpial_t');
else
    a=load(file_name);   
    wmpial_t=[a.wmpial_t(1:nsims); a.wmpial_t(length(a.wmpial_t)/2+1:length(a.wmpial_t)/2+nsims)];
end