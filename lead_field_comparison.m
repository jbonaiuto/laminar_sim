function lead_field_comparison(subj_info, varargin)

% Parse inputs
defaults = struct('surf_dir', 'D:/pred_coding/surf');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

% Pre-computed lead field matrices
load('C:\Users\jbonai\Dropbox\meg\layer_sim\lead_field_comparison\pialMU.mat');
load('C:\Users\jbonai\Dropbox\meg\layer_sim\lead_field_comparison\whiteMU.mat');

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

wm=gifti(white_mesh);
pial=gifti(pial_mesh);

pial_white_map=map_pial_to_white(white_mesh, pial_mesh, ...
        'mapType', 'link', 'origPial', orig_pial_mesh, ...
        'origWhite', orig_white_mesh);
    
white_pial_map=map_white_to_pial(white_mesh, pial_mesh, ...
    'mapType', 'link', 'origPial', orig_pial_mesh, ...
    'origWhite', orig_white_mesh);
    
pial_white_lf_dist=sqrt(sum((pialMU-whiteMU(pial_white_map,:)).^2,2));
white_pial_lf_dist=sqrt(sum((pialMU(white_pial_map,:)-whiteMU).^2,2));

pial_neighbor_lf_dist=zeros(1,size(pialMU,1));
for idx=1:size(pialMU,1)
    faces=union(union(find(pial.faces(:,1)==idx),find(pial.faces(:,2)==idx)), find(pial.faces(:,3)==idx));
    neighbor_faces=pial.faces(faces,:);
    neighbor_vertices=setdiff(neighbor_faces(:),idx);
    neighbor_lf_dist=mean(sqrt(sum((repmat(pialMU(idx,:),length(neighbor_vertices),1)-pialMU(neighbor_vertices,:)).^2,2)));
    pial_neighbor_lf_dist(idx)=neighbor_lf_dist;
end

white_neighbor_lf_dist=zeros(1,size(whiteMU,1));
for idx=1:size(whiteMU,1)
    faces=union(union(find(wm.faces(:,1)==idx),find(wm.faces(:,2)==idx)), find(wm.faces(:,3)==idx));
    neighbor_faces=wm.faces(faces,:);
    neighbor_vertices=setdiff(neighbor_faces(:),idx);
    neighbor_lf_dist=mean(sqrt(sum((repmat(whiteMU(idx,:),length(neighbor_vertices),1)-whiteMU(neighbor_vertices,:)).^2,2)));
    white_neighbor_lf_dist(idx)=neighbor_lf_dist;
end

pial_relative=pial_white_lf_dist./pial_neighbor_lf_dist';
white_relative=white_pial_lf_dist./white_neighbor_lf_dist';
disp(sprintf('pial=%.2f', length(find(pial_relative<1.0))./length(pial_white_lf_dist).*100));
disp(sprintf('white=%.2f', length(find(white_relative<1.0))./length(pial_white_lf_dist).*100));
lim=[min([pial_relative; white_relative]) prctile([pial_relative; white_relative],99.9)];
bins=linspace(lim(1),lim(2),100);
figure();
hold on;
[counts,centers]=hist(pial_relative,bins);
bar(centers,counts,'b');
[counts,centers]=hist(white_relative,bins);
bar(centers,counts,'r');
xlimit=xlim();
xlim([xlimit(1)/4 xlimit(2)]);
xlabel('Relative lead field difference RMS');
ylabel('Number of vertices');
legend('Pial','White matter');