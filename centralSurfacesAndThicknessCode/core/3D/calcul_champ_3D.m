function [V,T1,T2,B1,B2]=calcul_champ_3D(Vol,Eps_init,CA_GM,subiculum,opt)
%calcul_champ_3D Computes the 3D vector field as well as streamlines from inner to outer boundaries, streamlines from outer to inner boundaries and the lengths of streamlines starting/ending at any point x within the hippocampal domain
%
%   Arguments:
%   - Vol (structure): structure that defines the hippocampal domain (matrix and voxel size)
%   - Eps_init (matrix): defines the orientation of the normal vector to the ribbon boundaries
%   - CA_GM (structure): defines the domain of CA gray matter (matrix, voxel size)
%   - subiculum (structure): defines the domain of subiculum gray matter (matrix, voxel size)
%   - opt (structure): parameters of the thickness methods (anisotropy factors, kernel size)
%
%   Returns:
%   V (structure): vector field (matrices V.vx, V.vy, V.vz)
%   T1 (matrix): length at each x of Vol of the streamline from inner boundary to x
%   T2 (matrix): length at each x of Vol of the streamline from x to outer boundary
%   B1 (structure): all streamlines originating at any point within the domain and following
%       the vector field
%   B2 (structure): all streamlines originating at any point within the domain and following
%      the opposite vector field


% read anisotropy factors
ax=opt.aniso(1);
ay=opt.aniso(2);
az=opt.aniso(3);

% initialise vector field
%-- CA-SP
[uxCA,uyCA,uzCA]=initialise_champ_aniso(CA_GM,Eps_init.*CA_GM,opt.aniso);
%-- subiculum
[uxs,uys,uzs]=initialise_champ_aniso(subiculum,Eps_init.*subiculum,opt.aniso);
%-- CA-SP + subiculum
%[uxCA,uyCA,uzCA]=initialise_champ(CA_GM,Eps_init.*CA_GM);
%[uxs,uys,uzs]=initialise_champ(subiculum,Eps_init.*subiculum);
ux=uxs+uxCA;
uy=uys+uyCA;
uz=uzs+uzCA;
N=sqrt((ux/ax).^2+(uy/ay).^2+(uz/az).^2);
Z=find(Vol==1);
ux(Z)=ux(Z)./N(Z);
uy(Z)=uy(Z)./N(Z);
uz(Z)=uz(Z)./N(Z);
ux(isnan(ux))=0;
uy(isnan(uy))=0;
uz(isnan(uz))=0;

% compute vector field from domain, initial epsilon (direction of outer
% normal) and initialisation of the vector field
%[V]=calcule_V_frV3D(Vol,Eps_init,opt,ux,uy,uz);
[V]=calcule_V_frV3DNoyau(Vol,Eps_init,opt,ux,uy,uz);

% compute inner/outer streamlines for the vector field and the length of
% each streamline
%[B1,T1,B2,T2]=calcule_lignes3D(Vol,V.vx,V.vy,V.vz,opt.aniso);
[B1,T1,B2,T2]=calcule_lignes3D_aniso(Vol,V.vx,V.vy,V.vz,opt.aniso);

end
