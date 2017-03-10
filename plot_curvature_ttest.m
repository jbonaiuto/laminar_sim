function plot_curvature_ttest(subj_info, session_num, freq, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10, 'snr', 5, ...
    'surf_dir', 'd:\pred_coding\surf', 'mri_dir', 'd:\pred_coding\mri',...
    'curvature_type', 'k1', 'clip_vals', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
Nmeth=length(methodnames);
methods_to_plot=[1 4];

orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');

orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');

simmeshes={white_mesh,pial_mesh};
Nmesh=length(simmeshes);

pial=gifti(pial_mesh);
pial_inflated=gifti(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed_inflated.surf.gii'));
wm=gifti(white_mesh);
wm_inflated=gifti(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed_inflated.surf.gii'));

pial_curvature=compute_curvature(pial_mesh, 'curvature_type', params.curvature_type);
wm_curvature=compute_curvature(white_mesh, 'curvature_type', params.curvature_type);

[ax, metric_data]=plot_surface_metric(pial, pial_curvature, 'clip_vals', params.clip_vals, 'clabel', params.curvature_type, 'limits',[-.3 .3]);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);

[ax, metric_data]=plot_surface_metric(pial_inflated, pial_curvature, 'clip_vals', params.clip_vals, 'clabel', params.curvature_type, 'limits',[-.3 .3]);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);

[ax, metric_data]=plot_surface_metric(wm, wm_curvature, 'clip_vals', params.clip_vals, 'clabel', params.curvature_type, 'limits',[-.3 .3]);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);

[ax, metric_data]=plot_surface_metric(wm_inflated, wm_curvature, 'clip_vals', params.clip_vals, 'clabel', params.curvature_type, 'limits',[-.3 .3]);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);


nverts=size(wm.vertices,1);
rng(0);
simvertind=randperm(nverts); %% random list of vertex indices to simulate sources on

data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id, num2str(session_num), sprintf('f%d_%d_SNR%d_dipolemoment%d', freq(1), freq(2), params.snr, params.dipole_moment));

for i=1:length(methods_to_plot),       
    figure();
    methind=methods_to_plot(i);
    method=methodnames{methind};
    hold on;
    
    wmpial_t=get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
        white_mesh, orig_pial_mesh, orig_white_mesh);
    
    for simmeshind=1:Nmesh,    
        t_vals=wmpial_t((simmeshind-1)*params.nsims+1:simmeshind*params.nsims);
        if simmeshind==1
            t_vals=t_vals.*-1;
        end
        
        switch simmeshind
            case 1
                sim_curvature=wm_curvature(simvertind(1:params.nsims));
                color='r';
            case 2
                sim_curvature=pial_curvature(simvertind(1:params.nsims));
                color='b';
        end
        plot(sim_curvature,t_vals,['o' color]); 
    end
    hold off;
    xlabel(sprintf('Curvature (%s)',params.curvature_type))
    ylabel('t (correct-incorrect)');
    title(sprintf('t, %s',method)); 
end
