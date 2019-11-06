function create_median_surface(dim, Data_dir, Data_In_Prefix, Data_result_dir, Data_Out_Prefix, sigma) 
% create_median_surface Creates median surface with associated thickness map for either the 2D or the
% 3D version of the algorithm
%
% Arguments:
%   - dim (str): either '2D' or '3D' depending on the version of the algorithm the user wants to apply. Legacy, will always be 3D
%   - Data_dir: directory containing segmented images of hippocampal subregions in the BrainVisa IMA format. 
%        - the function expects each segmented subregion to be saved as an individual binary mask
%        - their respective filenames should contain a subregion identifier string: {'CA_GM','ca_gm','CA_SP','ca_sp','subiculum','SUBICULUM','CA_WM','ca_wm','sub_WM','SUB_WM','DG','dg','alveus','ALVEUS','fimbria','FIMBRIA'} . 
%   - Data_In_Prefix: legacy, not used
%   - Data_result_dir: output directory 
%   - Data_Out_Prefix (str): prefix of all the output data (images, meshes)
%   - outputPrefix: output file prefix
%   - sigma(int, optional): size of the kernel (influences the regularity of the vector field). Default value=10

sigma = str2num(sigma)

if strcmp(dim,'2D') == 1
  hippocamp_Seg2Skel(Data_dir,Data_result_dir,Data_Out_Prefix,sigma,1,1,1,1,cat(2,Data_Out_Prefix,'.txt'));
else
  hippocamp_Seg2Skel(Data_dir,Data_result_dir,Data_Out_Prefix,sigma,1,1,1,1,cat(2,Data_Out_Prefix,'.txt'));
end



