function plot_depth_free_energy(subj_info, session_num, freq, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'snr', 5, 'dipole_moment', 10, 'surf_dir', 'd:\pred_coding\surf',...
    'mri_dir', 'd:\pred_coding\mri', 'clip_vals', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
methods_to_plot=[1 4];
Nmeth=length(methodnames);

orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');

orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');

allmeshes=strvcat(white_mesh, pial_mesh);
Nmesh=size(allmeshes,1);

wm=gifti(allmeshes(1,:));
pial=gifti(allmeshes(2,:));
wm_inflated=gifti(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed_inflated.surf.gii'));
pial_inflated=gifti(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed_inflated.surf.gii'));

pial_white_map=map_pial_to_white(white_mesh, pial_mesh, ...
        'mapType', 'link', 'origPial', orig_pial_mesh, ...
        'origWhite', orig_white_mesh);

[pial_depth,HS]=compute_sulcal_depth(pial);
mapping=dsearchn(HS.vertices,wm.vertices);
wm_depth=sqrt(sum((wm.vertices-HS.vertices(mapping,:)).^2,2));

figure();
plot(pial_depth,wm_depth(pial_white_map),'.');
hold on
plot([0 60],[0 60],'r-');
xlabel('Pial depth');
ylabel('White matter depth');

[ax, metric_data]=plot_surface_metric(pial, pial_depth, 'clip_vals', params.clip_vals, 'limits', [0.5 50]);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);

[ax, metric_data]=plot_surface_metric(pial_inflated, pial_depth, 'clip_vals', params.clip_vals, 'limits', [0.5 50]);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);

[ax, metric_data]=plot_surface_metric(wm, wm_depth, 'clip_vals', params.clip_vals, 'limits', [0.5 50]);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);

[ax, metric_data]=plot_surface_metric(wm_inflated, wm_depth, 'clip_vals', params.clip_vals, 'limits', [0.5 50]);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);

nverts=size(wm.vertices,1);
rng(0);
simvertind=randperm(nverts); %% random list of vertex indices to simulate sources on

for i=1:length(methods_to_plot),       
    methind=methods_to_plot(i);
    figure();
    hold on;
    
    for simmeshind=1:Nmesh,    
        [path,file,ext]=fileparts(deblank(allmeshes(simmeshind,:)));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        
        data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),params.snr,params.dipole_moment));
        load(data_file);

        % F reconstructed on true - reconstructed on other
        % num simulations x number of folds
        truotherF=squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind)-allcrossF(simmeshind,1:params.nsims,(2-simmeshind)+1,methind));
        
        switch simmeshind
            case 1
                sim_depth=wm_depth(simvertind(1:params.nsims));
                color='r';
            case 2
                sim_depth=pial_depth(simvertind(1:params.nsims));
                color='b';
        end
        plot(sim_depth,truotherF,['o' color]);        
    end
    hold off;
    xlabel('Depth');
    ylabel('Free energy diff (correct-incorrect)');
    title(sprintf('Free energy, %s',methodnames{methind}));        
end
