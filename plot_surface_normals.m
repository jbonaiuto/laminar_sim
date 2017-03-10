function plot_surface_normals(surface)

norm    = spm_mesh_normals(struct('faces',surface.faces,'vertices',surface.vertices),true);  % normals to hipp surface
norm2   = norm*-1; %if you want them to point the other direction

%plot:
figure;
hold on;
quiver3(surface.vertices(1:10000,1), surface.vertices(1:10000,2),surface.vertices(1:10000,3),norm2(1:10000,1),norm2(1:10000,2),norm2(1:10000,3))
ax=gca;
plot_surface(surface,'ax',ax);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);

axis off;  set(gcf,'color','w'); hold on;%axis vis3d;
