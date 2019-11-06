function save_dim(file_dim,x_dim,y_dim,z_dim,dx,dy,dz,type,file_format)
dt=1;
bo='DCBA';
fid2 = fopen(file_dim,'w');
if fid2== -1, error(sprintf('[save_dim] Cannot open %s.',file_dim)); end
fprintf(fid2,'%d %d %d 1\n',x_dim,y_dim,z_dim);
fprintf(fid2,'-type %s\n',type);
fprintf(fid2,'-dx %f -dy %f -dz %f -dt %d\n',dx,dy,dz,dt);
fprintf(fid2,'-bo %s\n',bo);
fprintf(fid2,'-om %s\n',file_format);
fclose(fid2); 
if fid2 == -1, error(sprintf('[save_dim] Cannot close %s.',strcat(pathstr,name,'.dim'))); end