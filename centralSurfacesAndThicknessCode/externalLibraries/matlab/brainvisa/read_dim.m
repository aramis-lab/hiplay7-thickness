function [dim_mat,vox_size,type,file_format]=read_dim(file_dim)
fid2 = fopen(file_dim,'r');
if fid2== -1, error(sprintf('[load_ima] Cannot open %s.',file_dim)); end
%z_dim = fscanf(fid2,'%d',1)
x_dim = fscanf(fid2,'%d',1);
y_dim = fscanf(fid2,'%d',1);
%x_dim = fscanf(fid2,'%d',1)
z_dim = fscanf(fid2,'%d',1);
type=fscanf(fid2,'%s',1);
if isequal(type,'-type')==0
type=fscanf(fid2,'%s',1);
end
type=fscanf(fid2,'%s',1);
dx=fscanf(fid2,'%s',1);
dx=fscanf(fid2,'%f',1);
dy=fscanf(fid2,'%s',1);
dy=fscanf(fid2,'%f',1);
dz=fscanf(fid2,'%s',1);
dz=fscanf(fid2,'%f',1);
tmp1=fscanf(fid2,'%c',1);
tmp2=fscanf(fid2,'%c',1);
while ~strcmp([tmp1,tmp2],'om')
    tmp1=tmp2;
    tmp2=fscanf(fid2,'%c',1);
end
file_format=fscanf(fid2,'%s',1);
dim_mat=[x_dim,y_dim,z_dim];
vox_size=[dx,dy,dz];
fclose(fid2); 
if fid2 == -1, error(sprintf('[load_ima] Cannot close %s.',strcat(pathstr,name,'.dim'))); end