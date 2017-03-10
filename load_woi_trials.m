function woi_trials=load_woi_trials(woi_dir, prefix, woi, foi, nvertices)

% Load all pial data from woi
[files,dirs] = spm_select('List', fullfile(woi_dir, sprintf('%st%d_%d_f%d_%d*.gii', prefix, woi(1), woi(2), foi(1), foi(2))));
woi_trials=zeros(nvertices,size(files,1));
for t=1:size(files,1)
    filename=fullfile(woi_dir, sprintf('%st%d_%d_f%d_%d_1_%d.gii', prefix, woi(1), woi(2), foi(1), foi(2), t));
    trial_mesh=gifti(filename);
    woi_trials(:,t)=trial_mesh.cdata(:);
end