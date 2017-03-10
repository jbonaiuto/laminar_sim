function plot_paired_priors(subj_info, varargin)

% Parse inputs
defaults = struct('dir', 'c:/Users/jbonai/Dropbox/meg/layer_sim/paired_priors', 'surf_dir', 'd:\pred_coding\surf');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');

orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');
pialwhite_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed-ds_pial.hires.deformed.surf.gii');

pial=gifti(pial_mesh);
white=gifti(white_mesh);
pialwhite=gifti(pialwhite_mesh);

pial_white_map=map_pial_to_white(white_mesh, pial_mesh, ...
        'mapType', 'link', 'origPial', orig_pial_mesh, ...
        'origWhite', orig_white_mesh);
white_pial_map=map_white_to_pial(white_mesh, pial_mesh, ...
    'mapType', 'link', 'origPial', orig_pial_mesh, ...
    'origWhite', orig_white_mesh);
    
%% Setup simulation - number of sources, list of vertices to simulate on
nverts=size(white.vertices,1);
rng(0);
simvertind=randperm(nverts); %% random list of vertex indices to simulate sources on
Nsim=60; %% number of simulated sources on each surface

%% for MSP  or GS or ARD
% Number of patches as priors
% so use all vertices that will be simulated on, on white/pial surface
% (plus a few more) as MSP priors
%   wm sources         pial pairs for these sources               pial wources              wm pairs for these sources
Ip=[simvertind(1:Nsim) nverts+white_pial_map(simvertind(1:Nsim))' nverts+simvertind(1:Nsim) pial_white_map(simvertind(1:Nsim))'];

wm_pos=pialwhite.vertices([simvertind(1:Nsim) pial_white_map(simvertind(1:Nsim))'],:);
pial_pos=pialwhite.vertices([nverts+white_pial_map(simvertind(1:Nsim))' nverts+simvertind(1:Nsim)],:);

% Plot surface
ax=plot_surface(pialwhite, 'surface_alpha', 0.1);
fig=get(ax,'Parent');

% Plot coordinates
[x,y,z]=sphere();
rad=2.0;
for s=1:size(wm_pos,1)
    hp=surface(x.*rad+double(wm_pos(s,1)),y.*rad+double(wm_pos(s,2)),z.*rad+double(wm_pos(s,3)),...
       'FaceColor','r','EdgeColor','none','linestyle','none','FaceLighting','phong');
end
for s=1:size(pial_pos,1)
    hp=surface(x.*rad+double(pial_pos(s,1)),y.*rad+double(pial_pos(s,2)),z.*rad+double(pial_pos(s,3)),...
       'FaceColor','b','EdgeColor','none','linestyle','none','FaceLighting','phong');
end

set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-1.0 0.0 0.0]);
set(ax,'CameraPosition',[1760.839 -26.756 121.697]);

fig=get(ax,'Parent');
saveas(fig, fullfile(params.dir,'paired_prior_locations_side.png'), 'png');
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);
saveas(fig, fullfile(params.dir,'paired_prior_locations_topdown.png'), 'png');