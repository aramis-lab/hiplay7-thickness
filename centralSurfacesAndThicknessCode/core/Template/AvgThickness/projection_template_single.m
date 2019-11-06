function [EpProj]=projection_template_single(SkelSujPhi,EpSuj,Y,methode,sigma)
% [EpProj]=projection_template_single(SkelSujPhi,EpSuj,Y,methode,sigma)
% Project thickness values from a source to a target mesh
%
% Arguments
%   SkelSujPhi (cell): SkelSujPhi{1} is the source vertices (matrix)
%   EpSuj (cell): EpSuj{1} is the thickness values associated to the
%       source vertices (matrix)
%   Y (matrix): target vertices
%   methode (str): method used for the projection of thickness values.
%       Either 'gaussian' or 'wendland'
%   sigma (float): kernel size used when computing the projections
%
% Returns:
%   EpProj (matrix): thickness values associated to each vertex of the
%       target mesh (after projection of the thickness from source to
%       target)


N=numel(SkelSujPhi);

sigma2 = sigma^2;
EpProj=zeros(size(Y,1),N);

for i=1:N
    X = SkelSujPhi{i}';
    %X = X./repmat(aniso,size(X,1),1);
    E = EpSuj{i}';

    n = size(X,1);
   
    lambda = 10000;
    switch lower(methode)
        case{'gaussian'}

            D=(repmat(X(:,1),1,n)-repmat(X(:,1)',n,1)).^2+(repmat(X(:,2),1,n)-repmat(X(:,2)',n,1)).^2+...
                (repmat(X(:,3),1,n)-repmat(X(:,3)',n,1)).^2;

            K = exp(-D/sigma2)+(1/lambda)*eye(size(D));
            
            a = linsolve(K,E);
            m = size(Y,1);

            DXY=(repmat(Y(:,1),1,n)-repmat(X(:,1)',m,1)).^2+(repmat(Y(:,2),1,n)-repmat(X(:,2)',m,1)).^2 +...
                (repmat(Y(:,3),1,n)-repmat(X(:,3)',m,1)).^2;
            KXY = exp(-DXY/sigma2);
            %
        case{'wendland'}
            %
            D=sqrt((repmat(X(:,1),1,n)-repmat(X(:,1)',n,1)).^2+(repmat(X(:,2),1,n)-repmat(X(:,2)',n,1)).^2+...
                (repmat(X(:,3),1,n)-repmat(X(:,3)',n,1)).^2)/sigma;

            I=(D<=1);
            K=zeros(size(D));
            K(I) =(1-D(I)).^4.*(4*D(I)+1);
            K=K+(1/lambda)*eye(size(D));

            a = linsolve(K,E);
            m = size(Y,1);
            DXY=sqrt((repmat(Y(:,1),1,n)-repmat(X(:,1)',m,1)).^2+(repmat(Y(:,2),1,n)-repmat(X(:,2)',m,1)).^2 +...
                (repmat(Y(:,3),1,n)-repmat(X(:,3)',m,1)).^2)/sigma;


            KXY=zeros(size(DXY));
            IXY=(DXY<=1);
            KXY(IXY) =(1-DXY(IXY)).^4.*(4*DXY(IXY)+1);
        otherwise
            disp('Unknown method.')
    end


    EpProj(:,i) = KXY*a;


end
end








