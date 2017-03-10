function plot_curvature(surface_file, inflated, varargin)

% Parse inputs
defaults = struct('clip_vals', true, 'curvature_type', 'k1');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

curvature=compute_curvature(surface_file, 'curvature_type', params.curvature_type);

% Plot principal curvature k1
fig=figure('Position',[1 1 1200 400],'PaperUnits','points','PaperPosition',[1 1 600 200],'PaperPositionMode','manual');
ax=subplot(1,2,1);
plot_surface_metric(surface, curvature, 'ax', ax, 'clip_vals', params.clip_vals);
c=colorbar();
ylabel(c,params.curvature_type)
ax=subplot(1,2,2);
plot_surface_metric(inflated, curvature, 'ax', ax, 'clip_vals', params.clip_vals);
c=colorbar();
ylabel(c,params.curvature_type)