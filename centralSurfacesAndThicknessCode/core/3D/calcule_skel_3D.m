function [Skel,Col,SkelFace,SkelVert]=calcule_skel_3D(Vol,T1,T2,aniso)
% calcule_skel_3D Copmute central surface map (skeleton) as well as
% thickness map based on streamlines lengths
%
%   Arguments:
%   - Vol (structure): structure that defines the hippocampal domain (matrix and voxel size)
%   - T1 (matrix): length at each x of Vol of the streamline from inner boundary to x
%   - T2 (matrix): length at each x of Vol of the streamline from x to outer boundary
%
%   Returns:
%   - Skel (matrix): [skelVert_n, 3] matrix, central surface vertices of
%         the central surface mesh(mm Space)
%   - Col (vector): [skelVert_n] vector, thickness at each vertex
%   - SkelFace (matrix): [skelFace_n, 4] matrix, faces of the central
%         surface mesh
%   - SkelVert (matrix): [skelVert_n, 3] matrix, central surface vertices of
%         the central surface mesh(voxel Space)

fprintf('Compute central surface... ');

% Compute potential and thickness map
Mid=T1-T2;
Mid(~Vol)=NaN;
ThickTot=T1+T2;

% Generate mesh (coordinates of each voxel)
[l,m,n]=size(Vol);
ax=aniso(1);
ay=aniso(2);
az=aniso(3);
[A,B,C]=meshgrid(0:1/ay:(m-1)/ay,0:1/ax:(l-1)/ax,0:1/az:(n-1)/az);

% Compute isosurface          
[SkelFace,SkelVert]=isosurface(A,B,C,Mid,0);
Skel=SkelVert(:,[2 1 3]).*repmat(aniso,[size(SkelVert,1) 1])+1;

% interpolate thickness at each vertex
NS=size(Skel,1);
Col=zeros(NS,1);
for i=1:NS
    s=Skel(i,:);
    Col(i)=interp_tri(ThickTot,s(1),s(2),s(3));
end


end
