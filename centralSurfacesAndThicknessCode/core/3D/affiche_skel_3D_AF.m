function p=affiche_skel_3D_AF(SkelFace,SkelVert,Col,cmin,cmax,view_angle)
% function p=affiche_skel_3D_AF(SkelFace,SkelVert,Col,cmin,cmax,view_angle)
% 
% Display image of the skeleton mesh (defined by vertices, faces and
% thickness values at each vertex).
%
% Arguments:
%   SkelFace (matrix): faces of the skeleton mesh
%   SkelVert (matrix): vertices of the skeleton mesh
%   Col (matrix): thickness values associated with each vertex
%   cmin (float): minimum value for the colourmap
%   cmax (float): maximum value for the colourmap
%   view_angle (float): view angle for the displayed mesh in the Matlab
%       figure
%
% Returns:
%   p (patch object): patch object that contains the data for the
%       polygons in the Matlab figure

p=patch('Faces',SkelFace,'Vertices',SkelVert,'FaceVertexCData',Col,'EdgeColor','none','FaceColor','interp','CDataMapping','scaled');  
daspect([1 1 1])
view(view_angle) %[-90,90]
axis tight
caxis([cmin cmax])
camlight
lighting gouraud
set(gcf,'color','w','Renderer','zbuffer')

end
