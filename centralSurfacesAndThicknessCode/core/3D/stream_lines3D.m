function [B,Arret]=stream_lines3D(Vx,Vy,Vz,Start,Vol)
%stream_lines3D Computes the streamlines that follow a vector field within the domain for a range of starting positions
%
%   Arguments:
%   - Vx (matrix): x coordinates of the vector field
%   - Vy (matrix): y coordinates of the vector field
%   - Vz (matrix): z coordinates of the vector field
%   - Start (matrix): list of (x,y,z) starting positions for each streamline
%   - Vol (structure): structure that defines the hippocampal domain (matrix
%     and voxel size)
%
%   Returns:
%   - B (structure): all streamlines originating at any point within the domain
%     and following the vector field
%   - Arret (matrix): store for each streamline information about the
%     Arret(index,:)=[a,h,i,b] where
%     a: float difference between two consecutive points in the
%        integration scheme
%     h: step size
%     i: number of iterations used to reach the end of the integration
%     b: variable assessing whether the last point went outside the
%        domain or not


B=cell(1,size(Start,1));
Arret=zeros(size(Start,1),4);
n=size(Start,1);

for i=1:n
    % print current iteration to keep track of progress
    %i
    if mod(i,floor(n/100))==0
        disp(cat(2,int2str(i),'/',int2str(n)))
    
    end
    
    % compute the individual streamline corresponding to the current
    % starting position (Runge-Kutta integration)
    [S,arret]=RK_streamline3D(Vx,Vy,Vz,0.05,500,Start(i,:),Vol);
    
    % store result
    %-- streamline
    B{1,i}=S(:,[2, 1,3]);
    %-- ending position
    Arret(i,:)=arret;
    
end


end
