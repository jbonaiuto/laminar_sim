function plot_thickness(pial, wm, pial_inflated)

white_pial_map=map_white_to_pial(wm, pial);

dists=sum(sqrt((pial.vertices-wm.vertices(white_pial_map,:)).^2),2);
length(find(dists==0))

plot_surface_metric(pial_inflated, dists, 'clipped', true);