function [S,Arret] = RK_streamline3D(fx,fy,fz,h,N,P_init,A)
% RK_streamline3D solution of 1st order ODE using RK (Runge-Kutta) 2nd
% order with adaptive step size
%
%   Arguments:
%   - fx (matrix): 3D matrix. X-axis output of function f at the
%     position of each voxel within the matrix
%   - fy (matrix): 3D matrix. Y-axis output of function f at the
%     position of each voxel within the matrix
%   - fz (matrix): 3D matrix. Z-axis output of function f at the
%     position of each voxel within the matrix
%   - h (float): step size
%   - N (int): maximum number of iterations
%   - P_init (vector): 3-element vector starting position from which to
%     integrate the function f
%   - A (matrix): 3D matrix representing the domain within which
%     function f is defined
%     
%   Returns:
%   - output [t, y]  = solution vectors
%   - S (matrix): [N,3] matrix storing all the points that define the
%     streamline
%   - Arret (vector): 4-element vector to store information about the
%     parameters at the end of the Runge Kutta integration
%     Arret=[a,h,i,b] where
%     a: float difference between two consecutive points in the
%        integration scheme
%     h: step size
%     i: number of iterations used to reach the end of the integration
%     b: variable assessing whether the last point went outside the
%        domain or not



% initialise the streamline
S=zeros(N,3);
S(1,:)=P_init;


% define variables
%-- current iteration
i=1;
%-- float difference between two consecutive points in the integration
% scheme
a=1;
%-- threshold for the difference between two consecutive points
seuil=0.00001;
%-- variable to check whether the point found at the last iteration is
% still inside the domain
%b=A(round(S(1,1)),round(S(1,2)),round(S(1,3)));
%b=is_inside(A,S(1,1),S(1,2),S(1,3));
b=1;

%for i=1:N
% iterate while
% - max number of iterattions not reached
% - difference between consecutive points higher than threshold
% - streamline points still within domain
% - step size higher than threshold
while (((i<N)&(a>seuil))&(b==1)&(h>0.00001))

    % get current point
    pc=S(i,:);

    % successive interpolations
    k1x=interp_tri(fx,pc(1),pc(2),pc(3));
    k1y=interp_tri(fy,pc(1),pc(2),pc(3));
    k1z=interp_tri(fz,pc(1),pc(2),pc(3));

    k2x=interp_tri(fx,pc(1)+(h/2)*k1x,pc(2)+(h/2)*k1y,pc(3)+(h/2)*k1z);
    k2y=interp_tri(fy,pc(1)+(h/2)*k1x,pc(2)+(h/2)*k1y,pc(3)+(h/2)*k1z);
    k2z=interp_tri(fz,pc(1)+(h/2)*k1x,pc(2)+(h/2)*k1y,pc(3)+(h/2)*k1z);

    % define first point
    %S(i+1,:)=pc+ (h/2)*[k1x+k2x,k1y+k2y,k1z+k2z];
    X1=pc+ h*[k2x,k2y,k2z];
    
    % variable step size
    h2=h/2;


    % successive interpolations
    k1x=interp_tri(fx,pc(1),pc(2),pc(3));
    k1y=interp_tri(fy,pc(1),pc(2),pc(3));
    k1z=interp_tri(fz,pc(1),pc(2),pc(3));

    k2x=interp_tri(fx,pc(1)+(h2/2)*k1x,pc(2)+(h2/2)*k1y,pc(3)+(h2/2)*k1z);
    k2y=interp_tri(fy,pc(1)+(h2/2)*k1x,pc(2)+(h2/2)*k1y,pc(3)+(h2/2)*k1z);
    k2z=interp_tri(fz,pc(1)+(h2/2)*k1x,pc(2)+(h2/2)*k1y,pc(3)+(h2/2)*k1z);
    Xm=pc+ h2*[k2x,k2y,k2z];
    
    k1x=interp_tri(fx,Xm(1),Xm(2),Xm(3));
    k1y=interp_tri(fy,Xm(1),Xm(2),Xm(3));
    k1z=interp_tri(fz,Xm(1),Xm(2),Xm(3));

    k2x=interp_tri(fx,Xm(1)+(h2/2)*k1x,Xm(2)+(h2/2)*k1y,Xm(3)+(h2/2)*k1z);
    k2y=interp_tri(fy,Xm(1)+(h2/2)*k1x,Xm(2)+(h2/2)*k1y,Xm(3)+(h2/2)*k1z);
    k2z=interp_tri(fz,Xm(1)+(h2/2)*k1x,Xm(2)+(h2/2)*k1y,Xm(3)+(h2/2)*k1z);
    
    % define second point
    X2=Xm+ h2*[k2x,k2y,k2z];
    
    % adjust step size depending on how close first and second points
    % are
    delta=sqrt(sum((X2-X1).*(X2-X1)));
    
    if abs(delta) > 0.001
        h=h/2;
    else 
        if h*2<=0.2
            h=h*2;
        end
    end
    
    % append second point to the streamline
    %h
    % X2
    S(i+1,:)=X2;
    a=sum(abs(S(i+1,:)-S(i,:)));

    % check if point still inside the domain
    %b=1;
    %interp_tri(A,S(i+1,1),S(i+1,2),S(i+1,3))
     b=(interp_tri(A,S(i+1,1),S(i+1,2),S(i+1,3))>=0.5);
%     %if i>1
%         B=1-is_inside(A,S(i+1,1),S(i+1,2),S(i+1,3));
         if b==0
             % if outside domain, take previous point
             S(i+1,:)=Xm;
         end
    %end
     i = i+1;
    %b=b+B;
end
   
% process output
%-- remove empty points from the streamline
S(sum(S,2)==0,:)=[];
%-- store information about the end of the Runge-Kutta integration
Arret=[a,h,i,b];

end
