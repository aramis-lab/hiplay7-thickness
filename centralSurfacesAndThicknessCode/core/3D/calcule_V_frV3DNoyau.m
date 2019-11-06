function [V]=calcule_V_frV3DNoyau(Vol,Eps_init,Opt,ux,uy,uz)
% calcule_V_frV3DNoyau compute the transverse field that
% maximises the flow for a transverse field
%
%   Arguments:
%   - Vol (structure): structure that defines the hippocampal domain
%         (matrix and voxel size)
%   - Eps_init (matrix): initial orientation of the vector field
%         (inward/outward)
%   - Opt (structure): meta data and kernel size
%   - ux (matrix): approximate direction along X
%   - uy (matrix): approximate direction along Y
%   - uz (matrix): approximate direction along Z
%
%   Returns:
%   - V (structure): vector field, made of X, Y, Z components vx, vy
%         and vz. Also stores Vol, Eps (orientation of the vector field)
%         anisotropy (aniso), kernel size (sigma), S (variable
%         indicating how close successive vector fields are in the
%         transverse vector field iterative computation algorithm) and
%         niter (number of iterations taken to reach the output vector
%         field)

% default paramters
opt.sigma=5;
opt.aniso=[1 1 1];
if nargin>=3
    if (isfield(Opt,'sigma'))
        opt.sigma=Opt.sigma;
    end
    if (isfield(Opt,'niter'))
        opt.niter=Opt.niter;
    end
    if (isfield(Opt,'frV'))
        opt.frV=Opt.frV;
    end
    if (isfield(Opt,'aniso'))
        opt.aniso=Opt.aniso;
    end
    
end

% Iterative computation of the vector field
[vx vy vz Eps S h]=tvmflux_eg_frV_3D_anisoNoyau(Vol,ux,uy,uz,Eps_init,opt);


% Store all the output in structure
V.Vol=Vol;
V.Eps=Eps;
V.vx=vx;
V.vy=vy;
V.vz=vz;

V.sigma=opt.sigma;
V.niter=h;
V.S=S;
V.aniso=opt.aniso;


end
