function savemesh(filename,vertex,faces,normal,file_format)
%SAVEMESH Save a mesh in Brainvisa format
%  SAVEMESH(FILENAME,VERTEX,FACES,NORMAL) writes a mesh defined by
%  vertices VERTEX, triangles FACES and normals NORMAL in a Brainvisa
%  binary .mesh file FILENAME.
%  VERTEX and NORMAL  are (M x 3) arrays, while FACES is a (N x 3) array.
%  SAVEMESH(..,FILE_FORMAT) allows to specify the output file type: ascii
%  or binar. By default, FILE_FORMAT is 'binar'.
%
%  See also LOADMESH

% TODO % Handle .tri (ascii and binary)

%  Copyright (C) 2003 Guillaume Flandin
%  INRIA Sophia Antipolis / CEA-SHFJ
%
%  This program is free software; you can redistribute it and/or
%  modify it under the terms of the GNU General Public License
%  as published by the Free Software Foundation; either version 2
%  of the License, or any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation Inc, 59 Temple Pl. - Suite 330, Boston, MA 02111-1307, USA.

error(nargchk(4,5,nargin));

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
if fid == -1, error(sprintf('[savemesh] Cannot close %s.',filename)); end
