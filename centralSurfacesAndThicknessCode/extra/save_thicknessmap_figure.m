function save_thicknessmap_figure(thickmap, thickfig_folder, thickfig_prefix, colourbar_min, colourbar_max)
% Save the thickness map visualisation from two different views as a
% .png image.
% 
%
% Arguments:
%    thickmap (structure): structure containing the central surface volume
%    thickfig_folder (str): path to the folder where to output .png
%        visualisation
%    thickfig_prefix (str): prefix for the name of the files that
%        will be created inside thickfig_folder
%    colourbar_min (float): min value for the colourbar.
%    colourbar_max (float): max value for the colourbar.


% save view 1
thickfig_view1_path = strcat(thickfig_folder, '/', thickfig_prefix, '_view1.png') ;
figh = figure ;
p=affiche_skel_3D_AF(thickmap.SkelFace, thickmap.SkelVert, thickmap.SkelThick, colourbar_min, colourbar_max, [0, 90]); % add template thickness to the mean!
colorbar
colormap jet ;
cmap=colormap;
axis off; lighting gouraud;
saveas(figh, thickfig_view1_path) ;

% save view 2
thickfig_view2_path = strcat(thickfig_folder, '/', thickfig_prefix, '_view2.png') ;
figh = figure ;
p=affiche_skel_3D_AF(thickmap.SkelFace, thickmap.SkelVert, thickmap.SkelThick, colourbar_min, colourbar_max, [180, -90]); % add template thickness to the mean!
colorbar
colormap jet ;
cmap=colormap;
axis off; lighting gouraud;
saveas(figh, thickfig_view2_path) ;
end

