function L=longueur_ligne(A,dim,aniso)
% longueur_ligne Computes the length of a line
%
%   Arguments:
%   - A (cell): all streamlines
%   - dim (int): line dimension. 2 or 3.
%   - aniso (vector): [3] vector that defines the anisotropy factor
%         associated to the line
%
%   Returns:
%   - L (vector): vector containing all streamline lengths

% get number of streamlines
n=numel(A);


% initialisation
%-- initialise streamline lengths to l
L=zeros(n,1);
%-- check if anisotropy factor provided
if nargin==3
    aniso=aniso([2 1 3]);
else
    aniso=[1 1 1];
end

if dim==3
    for i=1:n
        % get current streamline and adjust with anisotropy factor
        ligne=A{1,i}./repmat(aniso,size(A{1,i},1),1);
        % check more than two points
        if numel(ligne)>2
            % sum distances between consecutive points
            l = sum(sqrt((ligne(2:end,1)-ligne(1:end-1,1)).^2+(ligne(2:end,2)-ligne(1:end-1,2)).^2+(ligne(2:end,3)-ligne(1:end-1,3)).^2));
            L(i)=l;
        else
            % less than two points
            L(i)=0;
        end
    end

else
    if dim==2
        for i=1:n
            % get current streamline and adjust with anisotropy factor
            ligne=A{1,i};
            if numel(ligne)>2
                % sum distances between consecutive points
                l = sum(sqrt((ligne(2:end,1)-ligne(1:end-1,1)).^2+(ligne(2:end,2)-ligne(1:end-1,2)).^2));
                L(i)=l;
            else
                % less than two points
                L(i)=0;
            end
        end
    end
end
end
