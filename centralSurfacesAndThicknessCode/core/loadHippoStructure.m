function [ structure ] = loadHippoStructure( dirData, listOfStructureNames, ignoreErrors )
% loadHippoStructure Loads an image file in brainvisa format '.ima' containing a given hippocampal structure
% 
% The file must be named according to one of the prefix values given in the
% argument listOfStructureNames. 
% If no file is found and ignoreErrors is 1 (default = 0), then the method raises an
% error. If loaded successfully, the function returns the image in the
% variable 'structure', containing the matrix (.mat). 
% If multiple files in the directory correspond to the searched filenames,
% only the first match will be considered.
%   
% Example : loadHippoStructure('/myhomedir/', {'*subiculum*', '*SUB*'})
% searches for files named '*subiculum*.ima' or '*SUB*.ima' within '/myhomedir/'
%
%   Arguments:
%   - dirData: directory containing segmented images of hippocampal subregions in the BrainVisa IMA format. 
%        - the function expects each segmented subregion to be saved as an individual binary mask
%        - their respective filenames should contain a subregion identifier string: {'CA_GM','ca_gm','CA_SP','ca_sp','subiculum','SUBICULUM','CA_WM','ca_wm','sub_WM','SUB_WM','DG','dg','alveus','ALVEUS','fimbria','FIMBRIA'} . 
%   - listOfStructureNames: list of potential names for the structures to be read (e.g., if looking for CA-SP, the
%        list would be '*CA_GM*', '*ca_gm*', '*CA_SP*', '*ca_sp*')
%   - ignoreErrors: Boolean. if true, will terminate code when a structure could not be read.
%        If false, will only display a warning when a structure could not be read.

    if ~exist('ignoreErrors')
        ignoreErrors = 1;
    end
    
    if dirData(end) ~= '/'
        dirData = strcat(dirData,'/');
    end

    % check if file exists for any of the potential names given in the
    % list of structure names
    structure=struct;
    files=[];
    nnames=numel(listOfStructureNames);
    name=1;
    while name <= nnames
        listOfStructureNames{name}
        strcat(dirData,listOfStructureNames{name},'.ima')
        files=cat(1,files,dir(strcat(dirData,listOfStructureNames{name},'.ima')));
        name=name+1;
        files = files(~cellfun('isempty', {files.date}))
        numel(files);
        if numel(files) > 0
            name=nnames+1;
        end
    end
    
    % if any file was found, read the first file that was found
    %files = files(~cellfun('isempty', {files.date}));
    if numel(files) > 0 
        temp=load_ima(strcat(dirData,files(1).name));
        %vox_size=temp.vox_size
        %dim_mat=temp.dim_mat
        %M=reshape(temp.mat,temp.dim_mat);
        
        % ensure the volume has binary values (0 and 1)
        structure.mat=temp.mat/max(temp.mat(:)); %M;
        structure.dim_mat=temp.dim_mat;
        
        structure.vox_size=temp.vox_size;
        structure.data_type=temp.data_type;
        
        structure
    else
        if ~ignoreErrors
            error(sprintf('No corresponding file found for the given hippocampal structure!'));
        else
            fprintf('No corresponding file found for the given hippocampal structure! Continuing anyway');
        end
    end
    
    return;
end

