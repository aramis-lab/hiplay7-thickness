function [ux,uy,uz]=initialise_champ(Vol,Eps_init)
%initialise_champ Initialise vector field from a binary domain and the
% initialisation of epsilon (direction of the outer normal at any boundary
% point).
% This function will work on isotropic input images.
%
%   Arguments:
%   - Vol (structure): structure that defines the hippocampal domain
%       (matrix and voxel size)
%   - Eps_init (matrix): defines the orientation of the normal vector to
%       the ribbon boundaries
%
%   Returns:
%   - ux (matrix): X coordinates of the vector field initialisation
%   - uy (matrix): Y coordinates of the vector field initialisation
%   - uz (matrix): Z coordinates of the vector field initialisation



% initialise the vector field from the gradient of distance to
% boundaries
%-- compute distance maps to boundaries
D1=bwdist(Eps_init==1);
D2=bwdist(Eps_init==-1);
%-- compute gradients for the distance maps
[fx,fy,fz]=gradient(D1);
[gx,gy,gz]=gradient(D2);
%-- compute vector field
uy=(fx.*(~(Eps_init==1))-gx.*(~(Eps_init==-1))).*Vol;
ux=(fy.*(~(Eps_init==1))-gy.*(~(Eps_init==-1))).*Vol;
uz=(fz.*(~(Eps_init==1))-gz.*(~(Eps_init==-1))).*Vol;
    
end
