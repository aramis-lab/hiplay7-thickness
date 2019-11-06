% plotmesh() - plot mesh defined by faces and vertex
%
% Usage: 
%     plotmesh(faces, vertex);
%
% Input:
%   faces   - array of N x 3. Each row defines a triangle. The 3 points
%             in each row are row indices in the matrix below.
%   vertex  - array of M x 3 points, (x = first colum; y=second colum
%             z=3rd column). Each row defines a point in 3-D.


% Copyright (C) May 6, 2003 Arnaud Delorme, SCCN/INC/UCSD,
% scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


function plotmesh(faces, vertex)
    figure; 

    FaceColor  = [.8 .55 .35]*1.1; % ~= ruddy Caucasian - pick your complexion!
    Lights = [-125  125  80; ...
              125  125  80; ...
              125 -125 125; ...
              -125 -125 125];    % default lights at four corners
    

    if any(any(faces == 0)), faces = faces+1; end;
    vertex(:,3) = -vertex(:,3);
    FCmap = [jet(64); FaceColor; FaceColor; FaceColor];
    colormap(FCmap)
    W = ones(1,size(vertex,1))*(size(FCmap,1)-1);
    p1 = patch('vertices', vertex, 'faces', faces, ...
               'FaceVertexCdata',W(:), 'FaceColor','interp');
    set(p1,'EdgeColor','none')
    

    for i = 1:size(Lights,1)
        hl(i) = light('Position',Lights(i,:),'Color',[1 1 1],...
                      'Style','infinite');
    end
    set(p1,'DiffuseStrength',.6,'SpecularStrength',0,...
           'AmbientStrength',.4,'SpecularExponent',5);
    axis equal
    view(18,8);
    rotate3d
    axis off;

