function [ ans ] = fileExists( filename )
%fileExists Checks if a given file exists
%   returns 1 in case of success, and 0 otherwise

    ans=0;

    fileindir = dir(filename);

    if numel(fileindir) > 0
        ans=1;
    end

end

