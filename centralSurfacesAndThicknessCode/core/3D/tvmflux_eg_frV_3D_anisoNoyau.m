function [vx vy vz Eps S h ]=tvmflux_eg_frV_3D_anisoNoyau(V,ux,uy,uz,Eps,opt)
% tvmflux_eg_frV_3D_anisoNoyau compute the transverse field that
% maximises the flow for a transverse field
%
%   Arguments:
%   - V (matrix): hippocampal domain
%   - ux (matrix): approximate direction along X
%   - uy (matrix): approximate direction along Y
%   - uz (matrix): approximate direction along Z
%   - Eps (matrix): initial orientation of the vector field
%         (inward/outward)
%   - opt (structure): meta data and kernel size
%
%   Returns:
%   - vx (matrix): X component of the vector field
%   - vy (matrix): Y component of the vector field
%   - vz (matrix): Z component of the vector field
%   - Eps (matrix): orientation of the vector field
%         (inward/outward)
%   - S (vector): records at each iteration whether new vector field is
%         close to previous vector field
%   - h (int): number of iterations taken to compute the vector field


% define default parameters in case those are not provided by the user
niter=100; % default max number of iterations
sigma=5; % default kernel size
aniso=[1 1 1]; % default anisotropy factor
frV=ones(size(V)); % 'frozen' voxels (where 0 values indicate freeze), for which no computations are done
if nargin==6
    if (isfield(opt,'sigma'))
        sigma=opt.sigma;
    end
    if (isfield(opt,'niter'))
        niter=opt.niter;
    end
    if (isfield(opt,'frV'))
        frV=opt.frV;
    end
    if (isfield(opt,'aniso'))
        aniso=opt.aniso;
    end
end

% define variables
%-- kernel-size relatd
l=ceil(sigma);
%%
%-- size of domain
[p,q,r]=size(V);
%-- anisotropy factors
ax=aniso(1);
ay=aniso(2);
az=aniso(3);
Det_a=ax*ay*az;


% adjust domain based on 'frozen' voxels
V(frV==0)=0;
% Compute normal vectors at the borders

% compute gradient of domain
%-- main computation
[By,Bx,Bz]=gradient(V);
%-- adjust with anisotropic factors
By=By.*V*ay/Det_a;
Bx=Bx.*V*ax/Det_a;
Bz=Bz.*V*az/Det_a;

%%

% initialise weights
%-- compute anisotropic kernel
filtre=filtre_aniso(ax,ay,az);
%-- deduce weights over the domain
Poids=(imfilter(V,filtre)).*V;
Poids(V~=0)=1./(Poids(V~=0));
PoidsC=sqrt(Poids);

% only keep gradient values inside the domain
Z=find(V~=0);
num=numel(Z);
Bx=Bx(Z); Bx=Bx(:);
By=By(Z); By=By(:);
Bz=Bz(Z); Bz=Bz(:);


e=Eps(Z);
e_init=e;
ind_e_init=find(e_init==0);

%% Initialise moments u
% Restrict to the domain
ux=ux(Z);
uy=uy(Z);
uz=uz(Z);
vx=zeros(p,q,r); vy=vx; vz=vx;
h=1;
s=1; % roughly: normalised sum of differences between current and previous vector field
seuil=1e-07; % Threshold

%%

while (h<=niter)||(s>seuil)
h
s
    % gradients multiplied by sign of outer normal
    eBx=e.*Bx;
    eBy=e.*By;
    eBz=e.*Bz;
    
    % keep current vector field in memory
    vxA=vx(Z); vyA=vy(Z); vzA=vz(Z);

    % compute new vector field over the domain
    vx(Z)=(ux/Det_a+eBx);
    vy(Z)=(uy/Det_a+eBy);
    vz(Z)=(uz/Det_a+eBz);
    

    % iterate until l (ceil(sigma)): update vector field based on
    % anisotropic weights
    for k=1:l
        vx=vx.*PoidsC;
        vy=vy.*PoidsC;
        vz=vz.*PoidsC;
        vx=imfilter(vx,filtre).*PoidsC;
        vy=imfilter(vy,filtre).*PoidsC;
        vz=imfilter(vz,filtre).*PoidsC;
    end
    
    % update ux, uy, uz
    Mux=vx(Z);
    Muy=vy(Z);
    Muz=vz(Z);

    normMu=sqrt((Mux/ax).^2+(Muy/ay).^2+(Muz/az).^2);
    
    ux=Mux./(normMu);
    uy=Muy./(normMu);
    uz=Muz./(normMu);
    
    ux(isnan(ux))=0;
    uy(isnan(uy))=0;
    uz(isnan(uz))=0;

    %ux=ux(:);uy=uy(:);uz=uz(:);

    ff=(rand(size(Mux))>0.5);

    %Me=(1-ff).*sign(Mux.*Bx+Muy.*By+Muz.*Bz)+ff.*e;
    %Me=sign(Mux.*Bx+Muy.*By+Muz.*Bz);
    Me=(1-ff).*sign(Mux.*Bx+Muy.*By+Muz.*Bz)+ff.*e;
    
    %numel(find((e_init~=Me)&(abs(e_init)==1)))

    e(ind_e_init)=Me(ind_e_init);
    %Eps(Z)=e;
    N=(cross([Mux,Muy,Muz],[vxA,vyA,vzA])).^2;%check if current vector field (Mu)
                                              %is close (colinear) to the vector
                                              %field at the previous stage (vA)
    s=sqrt(sum(N(:)))/num;
    S(h)=s;


    
    h=h+1;
end


%%

Eps=zeros(size(V));
Eps(Z)=e;

end

%%
function [filtre]=filtre_aniso(ax,ay,az)
% filtre_aniso compute the anisotropic kernel used in the computation
% of the vector field
%
%   Arguments:
%   - ax (float): anisotropy factor along X
%   - ay (float): anisotropy factor along Y
%   - az (float): anisotropy factor along Z
%
%   Returns:
%   - filtre (matrix): anisotropic kernel used in the computation of the
%     vector field
relativeErr = 1e-4;

% assign different edge weights based on the anisotropy factor
filtre=zeros(3,3,3);
if (abs(ax-1) < relativeErr*max(ax,1)) && (abs(ay-1) < relativeErr*max(ay,1)) && (az<1)
    L=1/az;
    
    alpha_x=L*L/(5*L*L+2);
    alpha_y=L*L/(5*L*L+2);
    alpha_z=1/(5*L*L+2);
    
    filtre([1,3],2,2)=alpha_x;
    filtre(2,[1,3],2)=alpha_y;
    filtre(2,2,[1,3])=alpha_z;
    
    filtre(2,2,2)=alpha_x;
    
else 
    if (abs(ax-1) < relativeErr*max(ax,1)) && (ay<1) && (abs(az-1) < relativeErr*max(az,1)) 
        L=1/ay;
        
        alpha_x=L*L/(5*L*L+2);
        alpha_y=1/(5*L*L+2);
        alpha_z=L*L/(5*L*L+2);
        
        filtre([1,3],2,2)=alpha_x;
        filtre(2,[1,3],2)=alpha_y;
        filtre(2,2,[1,3])=alpha_z;
    
        filtre(2,2,2)=alpha_x;
        
        
    else
        if (ax<1) && (abs(ay-1) < relativeErr*max(ay,1)) && (abs(az-1) < relativeErr*max(az,1))
            L=1/ax;
            
            alpha_x=1/(5*L*L+2);
            alpha_y=L*L/(5*L*L+2);
            alpha_z=L*L/(5*L*L+2);
            
            filtre([1,3],2,2)=alpha_x;
            filtre(2,[1,3],2)=alpha_y;
            filtre(2,2,[1,3])=alpha_z;
    
            filtre(2,2,2)=alpha_z;
            
        else 
            disp('incorrect anisotropy vector')
        end
    end
end


end
