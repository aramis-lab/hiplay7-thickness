function [ribbonVol_cropped,initEpsilon_cropped,CA_SP_cropped,SUB_cropped,metadata]=computeInnerOuterBoundaries2(segFilePath, outFilePath)
% computeInnerOuterBoundaries2 Read and compute volumes from the input data dir.
% All returned volumes are cropped around the hippocampal region of interest.
%
%   Arguments:
%   - segFilePath: directory containing segmented images of hippocampal subregions in the BrainVisa IMA format. 
%        - the function expects each each segmented subregion to be saved as an individual binary mask
%        - their respective filenames should contain a subregion identifier string: {'CA_GM','ca_gm','CA_SP','ca_sp','subiculum','SUBICULUM','CA_WM','ca_wm','sub_WM','SUB_WM','DG','dg','alveus','ALVEUS','fimbria','FIMBRIA'}. 
%   - outFilePath: output directory 
%
%   Returns:
%   - ribbonVol_cropped (structure): structure that defines the hippocampal domain (matrix and voxel size). Cropped around hippocampal ROI.
%   - initEpsilon_cropped (matrix): orientation of the vector normal to the hippocampal boundary
%   - CA_SP_cropped (structure): structure that defines the gray matter of CA (matrix and voxel size). Cropped around hippocampal ROI.
%   - SUB_cropped (structure): structure that defines the gray matter of subiculum (matrix and voxel size). Cropped around hippocampal ROI.
%   - metadata (structure): voxel size, volume dimensions, origin of the cropped ROI in the original (not cropped) ribbon volume

% tidy up paths
if segFilePath(end) ~= '/'
    segFilePath = strcat(segFilePath,'/');
end
if outFilePath(end) ~= '/'
    outFilePath = strcat(outFilePath,'/');
end


% Read volumes corresponding to hippocampal ribbon subregions 
%-- CA-SP
[CA_SP]=loadHippoStructure(segFilePath,{'*CA_GM*', '*ca_gm*', '*CA_SP*', '*ca_sp*'});
if ~isfield(CA_SP, 'mat')  
    return;
else
    save(cat(2,outFilePath,'CA_SP.mat'),'CA_SP');
end
%-- Subiculum
[SUB]=loadHippoStructure(segFilePath,{'*SUBICULUM*', '*Subiculum*', '*subiculum*'});
if ~isfield(SUB, 'mat')  
    return;
else
    save(cat(2,outFilePath,'SUB.mat'),'SUB');
end
%-- crash if dimensions are inconsistent
if ~(isequal(size(CA_SP.dim_mat),size(SUB.dim_mat)))
    return;
end

% compute hippocampal ribbon volume
ribbonVol=CA_SP.mat+SUB.mat;


% If inner and outer boundaries are available, load them 
[INNER]=loadHippoStructure(segFilePath,{'*bord_int*', '*inner_*'});
[OUTER]=loadHippoStructure(segFilePath,{'*bord_ext*', '*outer_*'});
if (isfield(INNER, 'mat') && isfield(OUTER, 'mat'))
    initEpsilon=sign(OUTER.mat-INNER.mat).*ribbonVol;
% Otherwise, read volumes corresponding to subregions in the inner or outer neighborhood of the ribbon and compute the boundaries
else
    % load Alveus
    [ALVEUS]=loadHippoStructure(segFilePath,{'*ALVEUS*', '*Alveus*', '*alveus*'});
    if ~isfield(ALVEUS, 'mat')  
        return;
    else
        save(cat(2,outFilePath,'ALVEUS.mat'),'ALVEUS');
    end

    % load CA_SRLM
    [CA_SRLM]=loadHippoStructure(segFilePath,{'*CA_WM*', '*Ca_wm*', '*ca_wm*', '*CA_SRLM*', '*Ca_srlm*', '*ca_srlm*'});
    if ~isfield(CA_SRLM, 'mat')  
        return;
    else
        save(cat(2,outFilePath,'CA_SRLM.mat'),'CA_SRLM');
    end

    % load DG
    [DG]=loadHippoStructure(segFilePath,{'*DG*', '*Dg*', '*dg*', '*GD*', '*Gd*', '*gd*'});
    if ~isfield(DG, 'mat')  
        return;
    else
        save(cat(2,outFilePath,'DG.mat'),'DG');
    end

    % load subiculum gray matter
    [SUB_WM]=loadHippoStructure(segFilePath,{'*SUBICULUM_WM*', '*Subiculum_WM*', '*Subiculum_wm*', '*subiculum_wm*', '*subiculum_WM*','*SUB_WM*', '*Sub_WM*', '*Sub_wm*', '*sub_wm*', '*sub_WM*'});
    if isfield(SUB_WM, 'mat')  
        save(cat(2,outFilePath,'SUB_WM.mat'),'SUB_WM');
    end

    % load fimbria
    [FIMBRIA]=loadHippoStructure(segFilePath,{'*FIMBRIA*', '*Fimbria*', '*fimbria*'});
    if isfield(FIMBRIA, 'mat') 
        save(cat(2,outFilePath,'FIMBRIA.mat'),'FIMBRIA');
    end

    % load uncus
    [UNCUS]=loadHippoStructure(segFilePath,{'*UNCUS*', '*Uncus*', '*uncus*'});
    if isfield(UNCUS, 'mat') 
        save(cat(2,outFilePath,'UNCUS.mat'),'UNCUS');
    end


    % check that subregions have compatible sizes
    if ~(isequal(size(CA_SP.dim_mat), size(CA_SRLM.dim_mat)) && isequal(size(CA_SP.dim_mat), size(ALVEUS.dim_mat)) && isequal(size(CA_SP.dim_mat), size(DG.dim_mat)))
        return;
    else
        if exist('SUB_WM.mat')
            if ~( isequal(size(CA_SP.dim_mat), size(SUB_WM.dim_mat)) )
                return;
            end
        end
    end


    % Compute inner and outer boundary
    filtre=conndef(3,'minimal');
    alveus_weight=imfilter(ALVEUS.mat,filtre);
    if exist('SUB_WM.mat')
        wm_weight=imfilter(CA_SRLM.mat+SUB_WM.mat,filtre);
    else
        wm_weight=imfilter(CA_SRLM.mat,filtre);
    end
    dg_weight=imfilter(DG.mat,filtre);

    initEpsilon=sign(alveus_weight-wm_weight).*ribbonVol-0.0001*(dg_weight>0).*ribbonVol;

    if exist('FIMBRIA.mat')
      if size(CA_SP.mat) == size(FIMBRIA.mat) 
          fimbria_weight=imfilter(FIMBRIA.mat,filtre);  
          initEpsilon=sign(initEpsilon+fimbria_weight.*ribbonVol);
      end
    end

    if exist('UNCUS.mat')
      if size(CA_SP.mat) == size(UNCUS.mat)
          hilum_weight=imfilter(UNCUS.mat,filtre);
          initEpsilon=sign(initEpsilon+hilum_weight.*ribbonVol);
      end
    end

end



Bord=bwperim(ribbonVol);

% /!\ Commented out as this is not necessary /!\
% Remaining unmarked subiculum points belong to the outer boundary
%Z=find((Bord==1)&(initEpsilon==0)&(SUB.mat==1)); 
%initEpsilon(Z)=1;


%% Initial \epsilon Eps_init
[maxval long_axis_dir]=max(CA_SP.vox_size);
last_slice=CA_SP.dim_mat(long_axis_dir);

i=last_slice;
while i > 0
    if long_axis_dir == 1
        if numel(find(ribbonVol(i,:,:) > 0)) > 0
	       initEpsilon(i,:,:)=0.0001; 
           last_slice=i;
           break; 
        end
    else
	    if long_axis_dir == 2
	        if numel(find(ribbonVol(:,i,:) > 0)) > 0
	            initEpsilon(:,i,:)=0.0001;
                last_slice=i;
                break; 
            end
	    else % long_axis_dir == 3
	        if numel(find(ribbonVol(:,:,i) > 0)) > 0
	            initEpsilon(:,:,i)=0.0001;
                last_slice=i;
                break;
            end
        end
    end
    i=i-1;
end

last_slice
long_axis_dir

if long_axis_dir == 1
    initEpsilon(last_slice,:,:)=0.0001;
else
    if long_axis_dir == 2
        initEpsilon(:,last_slice,:)=0.0001; %cancel extreme border (end of tail)
    else
        initEpsilon(:,:,last_slice)=0.0001;
    end
end

initEpsilon=initEpsilon.*ribbonVol;

%% Save boundaries, ribbon, and epsilon volumes
tempIma.dim_mat=CA_SP.dim_mat;
tempIma.vox_size=CA_SP.vox_size;
tempIma.mat=initEpsilon;
save_ima(strcat(outFilePath,'init_epsilon.ima'),tempIma,'DOUBLE','b');


tempIma.mat=ribbonVol;
save_ima(strcat(outFilePath,'hippo_ribbon.ima'),tempIma,'S16','b');


tempIma.mat=initEpsilon
tempIma.mat(find(initEpsilon > 0))=0;
tempIma.mat(find(initEpsilon < 0))=1;
min(min(min(tempIma.mat)))
save_ima(strcat(outFilePath,'inner_ribbon.ima'),tempIma,'S16','b');

tempIma.mat=initEpsilon
tempIma.mat(find(initEpsilon < 0))=0;
tempIma.mat(find(initEpsilon > 0.0001))=1;
save_ima(strcat(outFilePath,'outer_ribbon.ima'),tempIma,'S16','b');

% crop around ROI of the ribbon, to remove useless non-ribbon/non-boundary
% voxels

[ribbonVol_cropped.mat, initEpsilon_cropped.mat, CA_SP_cropped.mat, SUB_cropped.mat, origin]=cropRibbonVolumes(ribbonVol,initEpsilon,CA_SP.mat,SUB.mat);

metadata.dim_mat=CA_SP.dim_mat;
metadata.dim_mat_cropped=size(CA_SP_cropped.mat);
metadata.vox_size=CA_SP.vox_size;
metadata.origin=[0 0 0] + origin; %[minI-2,minJ-2,minK-2];
