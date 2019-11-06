function [ Vol_Ribbon,Vol_CA_SP,Vol_Sub,Vol_epsilon,metadata ] = generateInitialEpsilonFromRibbonBoundaries( dirData, dirResult, force )
% generateInitialEpsilonFromRibbonBoundaries Estimation of the initial epsilon values. If not yet available in dirData, the inner and outer
%                                  ribbon boundaries are computed based on the subregions segmentations, and saved in dirResult. 
%                                  The volume corresponding to a binary mask of the hippocampal ribbon is also saved in dirResult.
%   Arguments:
%   - dirData: directory containing segmented images of hippocampal subregions in the BrainVisa IMA format. 
%        - the function expects each each segmented subregion to be saved as an individual binary mask
%        - their respective filenames should contain a subregion identifier string: {'CA_GM','ca_gm','CA_SP','ca_sp','subiculum','SUBICULUM','CA_WM','ca_wm','sub_WM','SUB_WM','DG','dg','alveus','ALVEUS','fimbria','FIMBRIA'}. 
%   - dirResult: output directory 
%   - force: =1: to recompute all borders and initial epsilon; 
%            =0: to reuse previously computed ribbon borders and initial epsilon found under 'dirData' (if none of these is found, they are automatically computed)

% check input arguments
if nargin < 2
  error('input_example : dirData and dirResult are required inputs')
end
if nargin < 3
  force = 0
end
if dirData(end) ~= '/'
    dirData = strcat(dirData,'/');
end
if dirResult(end) ~= '/'
    dirResult = strcat(dirResult,'/');
end

% define names of the files of interest to be found in the input dir
%-- Hippocampal ribbon files
name_Roi = {'hippo_ribbon'};   
%-- Boundary files
name_inner = {'inner_ribbon'};   
name_outer = {'outer_ribbon'}; 
%-- initial epsilon file
name_eps = {'init_epsilon'};   
%-- hippocampal ribbon
file_ribbon = strcat(dirData,name_Roi{1},'.ima')
check_ribbon=fileExists(file_ribbon);
%-- initial epsilon (orientation of normal vector to the ribbon boundaries)
file_epsilon = strcat(dirData,name_eps{1},'.ima');
check_epsilon=fileExists(file_epsilon);



if force==1 || ~check_ribbon || ~(check_epsilon)
    % recompute
    sprintf('Computing boundaries and initial epsilon!')
    [Vol_Ribbon,Vol_epsilon,Vol_CA_SP,Vol_Sub,metadata]=computeInnerOuterBoundaries2(dirData, dirResult);
else
    % load from pre-existing data
    Vol_Ribbon=load_ima(file_ribbon);
    Vol_epsilon=load_ima(file_epsilon);
    Vol_CA_SP=loadHippoStructure(dirData,{'CA_GM', 'ca_gm', 'CA_SP', 'ca_sp'});
    Vol_Sub=loadHippoStructure(dirData,{'SUBICULUM', 'Subiculum', 'subiculum'});
    
    if (~isequal(size(Vol_Ribbon.dim_mat),size(Vol_epsilon.dim_mat))) || (~isequal(size(Vol_Ribbon.dim_mat),size(Vol_CA_SP.dim_mat))) || (~isequal(size(Vol_Ribbon.dim_mat),size(Vol_Sub.dim_mat)))
        error(sprintf('Existing hippocampal image sizes do not match!'));
        return;
    end

    metadata.vox_size=Vol_Ribbon.vox_size;
    metadata.dim_mat=Vol_Ribbon.dim_mat;

    % crop
    [Vol_Ribbon.mat, Vol_epsilon.mat, Vol_CA_SP.mat, Vol_Sub.mat, origin]=cropRibbonVolumes(Vol_Ribbon.mat,Vol_epsilon.mat,Vol_CA_SP.mat,Vol_Sub.mat);
    
    metadata.dim_mat_cropped=size(Vol_CA_SP.mat);
    metadata.origin=[0 0 0]+origin;
    
    metadata
end

end

