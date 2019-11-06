function [vertex, faces, normal, vertex_number, faces_number] = loadmesh(filename);
%LOADMESH Import Brainvisa mesh into Matlab
%  [VERTEX, FACES, NORMAL] = LOADMESH(FILENAME) loads a binary .mesh file FILENAME
%  and returns (M x 3) vertices in VERTEX, (M x 3) normals in NORMAL and (N x 3)
%  faces (triangles) in FACES.
%
%  See also SAVEMESH, PLOTMESH

% TODO % Handle .tri (ascii and binary)

%  Copyright (C) 2003 Denis Schwartz & Guillaume Flandin
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

error(nargchk(1,1,nargin));

fid = fopen(filename,'r');
if fid == -1, error(sprintf('[loadmesh] Cannot open %s.',filename)); end

[pathstr, name, format] = fileparts(filename);
if ~ismember(lower(format),{'.mesh'})
	error(sprintf('[loadmesh] Unknown format %s.',format));
end

[file_format, COUNT]       = fread(fid, 5, 'uchar');  %- 'ascii' or 'binar'

switch char(file_format)'
	case 'binar'
		[byte_swapping, COUNT]     = fread(fid, 1, 'uint32'); %- 'ABCD' or 'DCBA'
		ff = strcmp(dec2hex(byte_swapping),'41424344');
		if ~ff
			[fn, pm, mf] = fopen(1); %- machine format
			fclose(fid);
			if strmatch(mf,'ieee-le');
				fid = fopen(filename,'r','ieee-be');
			else
				fid = fopen(filename,'r','ieee-le');
			end
			[file_format, COUNT]   = fread(fid, 5, 'uchar');
			[byte_swapping, COUNT] = fread(fid, 1, 'uint32');
		end
		[arg_size, COUNT]          = fread(fid, 1, 'uint32'); %- length('VOID')
		[VOID, COUNT]              = fread(fid, arg_size, 'uchar'); %- VOID


		[polygon_dimension, COUNT] = fread(fid, 1, 'uint32'); %- 3 for triangles
		[mesh_time, COUNT]         = fread(fid, 1, 'uint32'); %- number of meshes

		vertex = cell(1,mesh_time);
		normals = cell(1,mesh_time);
		faces = cell(1,mesh_time);
		for i=1:mesh_time
			[mesh_step, COUNT]     = fread(fid, 1, 'uint32'); %- [0 ... mesh_time-1]

			%- Get vertices
			[vertex_number, COUNT] = fread(fid, 1, 'uint32');
			[vtx, COUNT]           = fread(fid, 3*vertex_number, 'float32');
			vertex{i} = reshape(vtx, 3, vertex_number)';

			%- Get normals
			[normal_number, COUNT] = fread(fid, 1, 'uint32');
			[nrml, COUNT]          = fread(fid, 3*normal_number, 'float32');
			normal{i} = reshape(nrml, 3, normal_number)';

			[arg_size, COUNT]      = fread(fid, 1, 'uint32'); %- no data ('VOID')

			%- Get faces
			[faces_number, COUNT]  = fread(fid, 1, 'uint32');
			[fcs, COUNT] = fread(fid, polygon_dimension*faces_number, 'uint32');
			faces{i} = reshape(fcs, polygon_dimension, faces_number)';
		end
	case 'ascii'
		VOID = fscanf(fid,'%s',1);
		polygon_dimension = fscanf(fid,'%d',1);
		mesh_time = fscanf(fid,'%d',1);
		
		for i=1:mesh_time
			mesh_step = fscanf(fid,'\n%d',1);
			
			vertex_number = fscanf(fid,'\n%d\n',1);
			vtx = fscanf(fid,'(%f ,%f ,%f) ',3*vertex_number);
			vertex{i} = reshape(vtx, 3, vertex_number)';
			
			normal_number = fscanf(fid,'\n%d\n',1);
			nrml = fscanf(fid,'(%f ,%f ,%f) ',3*normal_number);
			normal{i} = reshape(nrml, 3, normal_number)';
			
			arg_size = fscanf(fid,'\n%d\n',1);
			
			faces_number = fscanf(fid,'\n%d\n',1);
			fcs = fscanf(fid,'(%d ,%d ,%d) ',polygon_dimension*faces_number);
			faces{i} = reshape(fcs, polygon_dimension, faces_number)';
		end
end

if mesh_time == 1
	vertex = vertex{1};
	normal = normal{1};
	faces = faces{1};
end

fclose(fid); 
if fid == -1, error(sprintf('[loadmesh] Cannot close %s.',filename)); end
