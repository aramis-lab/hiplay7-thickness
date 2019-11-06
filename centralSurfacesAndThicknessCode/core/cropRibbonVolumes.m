function [ribbonVol_cropped, initEpsilon_cropped, CA_SP_cropped, SUB_cropped, origin] = cropRibbonVolumes(ribbonVol,initEpsilon,CA_SP,SUB)
% cropRibbonVolumes Crops image matrices tightly around voxels of the hippocampal ribbon and its
% boundaries and adds a 1-voxel border to the final images
%   Arguments:
%   - ribbonVol (structure): structure that defines the hippocampal domain (matrix and voxel size).
%   - initEpsilon_cropped (matrix): orientation of the vector normal to the hippocampal boundary
%   - CA_SP_cropped (structure): structure that defines the gray matter of CA (matrix and voxel size).
%   - SUB_cropped (structure): structure that defines the gray matter of subiculum (matrix and voxel size).
%
%   Returns:
%   - ribbonVol_cropped (structure): structure that defines the hippocampal domain (matrix and voxel size). Cropped around hippocampal ROI.
%   - initEpsilon_cropped (matrix): orientation of the vector normal to the hippocampal boundary
%   - CA_SP_cropped (structure): structure that defines the gray matter of CA (matrix and voxel size). Cropped around hippocampal ROI.
%   - SUB_cropped (structure): structure that defines the gray matter of subiculum (matrix and voxel size). Cropped around hippocampal ROI.
%   - origin (vector): 3-element vector, starting position in the original (non-cropped) domain of the cropped domain

% find crop boundaries from ribbon volumes and init Epsilon
mininiteps=min(min(min(initEpsilon)))
maxiniteps=max(max(max(initEpsilon)))

Vol_sum = ribbonVol + initEpsilon;
Z=find(Vol_sum~=0); % indices of non-zero voxels
[I,J,K]=ind2sub(size(Vol_sum),Z); % (x,y,z) indices of non-zero voxels

minI = min(I);
minJ = min(J);
minK = min(K);

maxI = max(I);
maxJ = max(J);
maxK = max(K);

% crop all volumes using crop boundaries
%-- initialse crops
ribbonVol_cropped = zeros(maxI+3-minI,maxJ+3-minJ,maxK+3-minK);
initEpsilon_cropped = zeros(maxI+3-minI,maxJ+3-minJ,maxK+3-minK); 
CA_SP_cropped = zeros(maxI+3-minI,maxJ+3-minJ,maxK+3-minK);
SUB_cropped = zeros(maxI+3-minI,maxJ+3-minJ,maxK+3-minK);

size(CA_SP)
size(CA_SP_cropped)

%-- fill crops
ribbonVol_cropped(2:end-1,2:end-1,2:end-1) = ribbonVol(minI:maxI,minJ:maxJ,minK:maxK);
initEpsilon_cropped(2:end-1,2:end-1,2:end-1) = initEpsilon(minI:maxI,minJ:maxJ,minK:maxK);
CA_SP_cropped(2:end-1,2:end-1,2:end-1) = CA_SP(minI:maxI,minJ:maxJ,minK:maxK);
SUB_cropped(2:end-1,2:end-1,2:end-1) = SUB(minI:maxI,minJ:maxJ,minK:maxK);

% define origin point
%origin=[minI-2,minJ-2,minK-2];
origin=[minI-1,minJ-1,minK-1];

mininiteps=min(min(min(initEpsilon)))
maxiniteps=max(max(max(initEpsilon)))

end

