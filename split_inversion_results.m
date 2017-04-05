function split_inversion_results(data_dir, varargin)

% Parse inputs
defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

[files,dirs] = spm_select('List', data_dir, '^*\.gii');
for f=1:size(files,1)
    filename=deblank(files(f,:));
    [path prefix ext]=fileparts(filename);
    trial_mesh=gifti(fullfile(data_dir,filename));
    n_combined_vertices=length(trial_mesh.cdata);
    n_vertices=round(n_combined_vertices/2);

    write_metric_gifti(fullfile(data_dir,['pial_' prefix]), trial_mesh.cdata(n_vertices+1:end,:));

    write_metric_gifti(fullfile(data_dir,['white_' prefix]), trial_mesh.cdata(1:n_vertices,:));

    delete(fullfile(data_dir,[prefix '.dat']),fullfile(data_dir,[prefix '.gii']));
end
