function [B1,Thick1,B2,Thick2]=calcule_lignes3D_aniso(Vol,Vx,Vy,Vz,aniso)
%calcule_lignes_3D_aniso Computes the streamlines from inner to outer boundaries, streamlines from outer to inner boundaries and the lengths of streamlines starting/ending at any point x within the hippocampal domain
%
%   Arguments:
%   - Vol (structure): structure that defines the hippocampal domain (matrix and voxel size)
%   - Vx (matrix): x coordinates of the vector field
%   - Vy (matrix): y coordinates of the vector field
%   - Vz (matrix): z coordinates of the vector field
%   - aniso (vector): anisotropy factor
%
%   Returns:
%   B1 (structure): all streamlines originating at any point within the domain and following
%       the vector field
%   Thick1 (matrix): length at each x of Vol of the streamline from inner boundary to x
%   B2 (structure): all streamlines originating at any point within the domain and following
%      the opposite vector field
%   Thick2 (matrix): length at each x of Vol of the streamline from x to outer boundary


% define starting positions as all the points within the volume
Z=find(Vol);
[I,J,K]=ind2sub(size(Vol),Z);
Start=[I,J,K];

% apply anisotropy factor to the vector field
ax=aniso(1);
ay=aniso(2);
az=aniso(3);
normV=sqrt((Vx/ax).^2+(Vy/ay).^2+(Vz/az).^2);
Vx=Vx./normV.*Vol;
Vy=Vy./normV.*Vol;
Vz=Vz./normV.*Vol;
Vx(isnan(Vx))=0;
Vy(isnan(Vy))=0;
Vz(isnan(Vz))=0;


% get B1: integrate along the vector field from each starting point
[B1,arret1]=stream_lines3D(Vx,Vy,Vz,Start,Vol);

% computing corresponding thicknesses
L1=longueur_ligne(B1,3,aniso);
Thick1=zeros(size(Vol));
Thick1(Z)=L1;



% get B2: integrate along minus the vector field from each starting point
Vxm=-Vx;Vym=-Vy;Vzm=-Vz;
[B2,arret2]=stream_lines3D(Vxm,Vym,Vzm,Start,Vol);

% computing corresponding thicknesses
L2=longueur_ligne(B2,3,aniso);
Thick2=zeros(size(Vol));
Thick2(Z)=L2;


end
