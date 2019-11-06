function savetex(filename,texture)
%SAVEMESH Save a texture in Brainvisa format



[pathstr, name, ext] = fileparts(filename);
if ~ismember(lower(ext),{'.tex'})
	error(sprintf('[savetex] Unknown format %s.',ext));
end

fid = fopen(filename,'w');
if fid == -1, error(sprintf('[savemesh] Cannot open %s.',filename)); end

nb_points=length(texture);
if size(texture,2)==1
    texture=texture';
end

fprintf(fid, 'ascii\n');
fprintf(fid, 'FLOAT\n');

fprintf(fid, '1\n');
fprintf(fid, '0\n');
fprintf(fid, '%d\n', nb_points);

fprintf(fid, '%f ', texture);
    


fclose(fid); 
if fid == -1, error(sprintf('[savemesh] Cannot close %s.',filename)); end
