function [ux,uy,uz]=initialise_champ_aniso(Vol,Eps_init,aniso)
%initialise_champ_aniso Initialise vector field from a binary domain,
% and the initialisation of epsilon (direction of the outer normal at
% any boundary point).
% This function has been designed to take the anisotropy factor of the
% input image into account.
%
%   Arguments:
%   - Vol (structure): structure that defines the hippocampal domain
%       (matrix and voxel size)
%   - Eps_init (matrix): defines the orientation of the normal vector to
%       the ribbon boundaries
%   - aniso (vector): 3-element vector defining the anisotropy factor
%       of the input image
%
%   Returns:
%   - ux (matrix): X coordinates of the vector field initialisation
%   - uy (matrix): Y coordinates of the vector field initialisation
%   - uz (matrix): Z coordinates of the vector field initialisation

% read anisotropy factor
ax=aniso(1);
ay=aniso(2);
az=aniso(3);

if sum([ax ay az]==[1 1 1])~=3
    % image is anisotropic

    [l,m,n]=size(Vol);

    % compute grids
    [X,Y,Z]=meshgrid(1:ay:m,1:ax:l,1:az:n);
    %-- define 'anisotropic' grid
    [x,y,z]=meshgrid(1:1/ay:m/ay,1:1/ax:l/ax,1:1/az:n/az);

    Eps_init2=interp3(Eps_init,X,Y,Z,'nearest');

    % initialise the vector field from the gradient of distance to
    % boundaries
    %-- compute distance maps to boundaries
    D1=bwdist(Eps_init2==1);
    D2=bwdist(Eps_init2==-1);
    %-- compute gradients for the distance maps
    [Fx,Fy,Fz]=gradient(D1);
    [Gx,Gy,Gz]=gradient(D2);
    %-- interpolate in 'anisotropic' grid
    fx=interp3(Fx,x,y,z,'nearest');
    fy=interp3(Fy,x,y,z,'nearest');
    fz=interp3(Fz,x,y,z,'nearest');
    gx=interp3(Gx,x,y,z,'nearest');
    gy=interp3(Gy,x,y,z,'nearest');
    gz=interp3(Gz,x,y,z,'nearest');
    %-- compute vector field
    uy=(fx.*(~(Eps_init==1))-gx.*(~(Eps_init==-1))).*Vol;
    ux=(fy.*(~(Eps_init==1))-gy.*(~(Eps_init==-1))).*Vol;
    uz=(fz.*(~(Eps_init==1))-gz.*(~(Eps_init==-1))).*Vol;

else 
    % image is isotropic. Use isotropic version of code
    [ux,uy,uz]=initialise_champ(Vol,Eps_init);
end
    
    
end
