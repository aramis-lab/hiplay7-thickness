function In=interp_tri(f,x,y,z)
% interp_tri trilinear interpolation
%
%   Arguments:
%   - f (matrix): values of the function to interpolate over a 3D grid
%   - x (float): X coordinate of the point at which to interpolate the
%     function
%   - y (float): Y coordinate of the point at which to interpolate the
%     function
%   - z (float): Z coordinate of the point at which to interpolate the
%     function
%
%   Returns:
%   - In (float): interpolated value of f at (x,y,z)

% define float difference to bottom-left-south voxel
xd=x-floor(x);
yd=y-floor(y);
zd=z-floor(z);

% define bottom-left-south voxel
xi=floor(x);
yi=floor(y);
zi=floor(z);

% i1=f(xi,yi,zi)*(1-zd)+f(xi,yi,zi+1)*zd;
% i2=f(xi,yi+1,zi)*(1-zd)+f(xi,yi+1,zi+1)*zd;
% j1=f(xi+1,yi,zi)*(1-zd)+f(xi+1,yi,zi+1)*zd;
% j2=f(xi+1,yi+1,zi)*(1-zd)+f(xi+1,yi+1,zi+1)*zd;

% linear interpolation
I=f(xi:xi+1,yi:yi+1,zi)*(1-zd)+f(xi:xi+1,yi:yi+1,zi+1)*zd;

%W=I*[1-yd;yd];
W=double(I)*[1-yd;yd];

In=W'*[1-xd;xd];

end


