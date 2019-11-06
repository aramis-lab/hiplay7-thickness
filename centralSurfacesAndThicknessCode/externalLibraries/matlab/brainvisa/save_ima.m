function save_ima(filename,ima,type,file_format)

fid = fopen(filename,'w');
if fid == -1, error(sprintf('[write_ima] Cannot open %s.',filename)); end

[pathstr, name, format] = fileparts(filename)
file_dim=strcat(pathstr,'/',name,'.dim');
save_dim(file_dim,ima.dim_mat(1),ima.dim_mat(2),ima.dim_mat(3),ima.vox_size(1),ima.vox_size(2),ima.vox_size(3),type,file_format);
switch type
    case 'S16'
        switch file_format
            case {'binar', 'b'}
                fwrite(fid,ima.mat, 'uint16');
            case {'ascii', 'a'}
                fprintf('rien ne se passe!\n');
            otherwise
                error('[to_ima] Unknown output format.');
        end
    case 'FLOAT'
        switch file_format
            case {'binar', 'b'}
                fwrite(fid,ima.mat, 'float');
            case {'ascii', 'a'}
                fprintf('rien ne se passe!\n');
            otherwise
                error('[to_ima] Unknown output format.');
        end
    case 'DOUBLE'
        switch file_format
            case {'binar', 'b'}
                fwrite(fid,ima.mat, 'double');
            case {'ascii', 'a'}
                fprintf('rien ne se passe!\n');
            otherwise
                error('[to_ima] Unknown output format.');
        end     
end

fclose(fid);
if fid == -1, error(sprintf('[to_ima] Cannot close %s.',filename)); end
