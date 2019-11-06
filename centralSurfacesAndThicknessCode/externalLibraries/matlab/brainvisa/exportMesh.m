function exportMesh( filename,vertex,faces,texture,normal,file_format)
error(nargchk(3,6,nargin));

if (~exist('normal','var'))
    normal = vertex;
end

if (~exist('file_format','var'))
    file_format = 'binar';
end

faces = faces - ones(size(faces));

[pathstr, name, ext] = fileparts(filename);
if ~ismember(lower(ext),{'.mesh'})
	error(sprintf('[savemesh] Unknown format %s.',ext));
end

fid = fopen(filename,'w');
if fid == -1, error(sprintf('[savemesh] Cannot open %s.',filename)); end

if nargin == 4
	file_format = 'binar';
else
	file_format = lower(file_format);
end

if iscell(vertex)
	nb_mesh = length(vertex);
else
	nb_mesh = 1;
	vertex = { vertex };
	faces  = { faces  };
	normal = { normal };
end

switch file_format
	case {'binar', 'b'}
		fwrite(fid, 'binar', 'uchar');
		fwrite(fid, hex2dec('41424344'), 'uint32');
		fwrite(fid, 4, 'uint32');
		fwrite(fid, 'VOID', 'uchar');

		fwrite(fid, size(faces{1},2), 'uint32');
		fwrite(fid, nb_mesh, 'uint32');

		for i=1:nb_mesh
			fwrite(fid, i-1, 'uint32');
	
			fwrite(fid, size(vertex{i},1), 'uint32');
			fwrite(fid, vertex{i}', 'float32');
	
			fwrite(fid, size(normal{i},1), 'uint32');
			fwrite(fid, normal{i}', 'float32');
	
			fwrite(fid, 0, 'uint32'); %- no data ('VOID')
	
			fwrite(fid, size(faces{i},1), 'uint32');
			fwrite(fid, faces{i}', 'uint32');
		end
	case {'ascii', 'a'}
		fprintf(fid, 'ascii\n');
		fprintf(fid, 'VOID\n');
		
		fprintf(fid, '%d\n', size(faces{1},2));
		fprintf(fid, '%d\n', nb_mesh);
		
		for i=1:nb_mesh
			fprintf(fid, '%d\n', i-1);
			
			fprintf(fid, '\n%d\n', size(vertex{i},1));
			fprintf(fid, '(%f ,%f , %f) ', vertex{i}');
			
			fprintf(fid, '\n%d\n', size(normal{i},1));
			fprintf(fid, '(%f ,%f ,%f) ', normal{i}');
			
			fprintf(fid, '\n0'); %- no data ('VOID')
			
			fprintf(fid, '\n%d\n', size(faces{i},1));
			fprintf(fid, '(%d ,%d ,%d) ', faces{i}');
		end
	otherwise
		error('[savemesh] Unknown output format.');
end

fclose(fid); 
if fid == -1, error(sprintf('[savemesh] Cannot close %s.',filename));end

if (exist('texture','var'))

    filename = strcat(pathstr,'/',name,'.tex');

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
end 

end

