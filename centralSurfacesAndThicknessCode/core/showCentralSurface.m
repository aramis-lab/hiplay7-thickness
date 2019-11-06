function showCentralSurface(matFilePath, outputFilePath, prefix)
% Given a path to a .mat filename with results from a central surface and thickness computation, display in matlab the thickness map plotted over the surface.
% Snapshots of the views are saved in the directory given by outputFilePath with files prefixed as "prefix"
% Ex.: showCentralSurface('./mymatfile.mat', './results/', 'central_surface_image')
results=load(matFilePath)

if ~exist('results', 'var')
  return
end

minv=min(results.SkelThick)
maxv=max(results.SkelThick)

figh=figure
p=affiche_skel_3D_AF(results.SkelFace,results.SkelVert,results.SkelThick, minv, maxv, [90,-90]);
axis off;
saveas(figh, strcat(outputFilePath, prefix, 'view1.png'))

p=affiche_skel_3D_AF(results.SkelFace,results.SkelVert,results.SkelThick, minv, maxv, [-90,90]);
colorbar;
axis off;
saveas(figh, strcat(outputFilePath, prefix, 'view2.png'))
