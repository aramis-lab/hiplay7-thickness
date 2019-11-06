function create_anisotropic_segmentations(input_segmentations_path, output_segmentations_path, subsample_factor)
% create_anisotropic_segmentations Starting from quasi-isotropic
% segmentations, generate 'simulated' segmentations from anisotropic MR
% data. Instead of subsampling the quasi-isotropic MR data to
% anisotropic MR data and carrying out segmentation of the anisotropic
% data, this subsamples directly the segmentations from the
% quasi-isotropic data.
%
% Arguments:
%   - input_segmentations_path (str): path to the folder that contains
%         the segmentation of the quasi-isotropic data. Should contain
%         the following files: [ALVEUS].dim/ima, [CA_GM].dim/ima,
%         [CA_WM].dim/ima, [SUBICULUM].dim/ima, [SUB_WM].dim/ima,
%         [DG].dim/ima
%         Where [XYZ] indicates any reasonable variations of 'xyz'
%         (e.g., all capital letters, all small letters, mixed
%         caps/small, etc.)
%   - output_segmentations_path (str): path to the folder that contains
%         the 'simulated' segmentations from anisotropic data. After
%         the function run, will contain the following files:
%         alveus.dim/ima, ca_gm.dim/ima, ca_wm.dim/ima
%         subiculum.dim/ima, sub_wm.dim/ima, dg.dim/ima
%   - subsample_factor (int): subsampling factor f, where the
%         anisotropic 'simulated' segmentations are produced by taking
%         one slice out of 'f'. This, along the anisotropy factor of the
%         original quasi-isotropic segmentation, defines the anisotropy
%         factor of the simulated segmentation.

% convert strings to number when required
subsample_factor = str2num(subsample_factor) ;

% load input segmentations
disp(' ') ;
disp('load input segmentations: alveus, CA_GM, CA_WM, subiculum, SUB_WM, DG') ;
[in_seg.alveus] = loadHippoStructure(input_segmentations_path, {'*ALVEUS*', '*Alveus*', '*alveus*'}) ;
[in_seg.ca_gm] = loadHippoStructure(input_segmentations_path, {'*CA_GM*', '*ca_gm*', '*CA_SP*', '*ca_sp*'}) ;
[in_seg.ca_wm] = loadHippoStructure(input_segmentations_path, {'*CA_WM*', '*Ca_wm*', '*ca_wm*', '*CA_SRLM*', '*Ca_srlm*', '*ca_srlm*'}) ;
[in_seg.subiculum] = loadHippoStructure(input_segmentations_path, {'*SUBICULUM*', '*Subiculum*', '*subiculum*'}) ;
[in_seg.sub_wm] = loadHippoStructure(input_segmentations_path, {'*SUBICULUM_WM*', '*Subiculum_WM*', '*Subiculum_wm*', '*subiculum_wm*', '*subiculum_WM*','*SUB_WM*', '*Sub_WM*', '*Sub_wm*', '*sub_wm*', '*sub_WM*'}) ;
[in_seg.dg] = loadHippoStructure(input_segmentations_path, {'*DG*', '*Dg*', '*dg*', '*GD*', '*Gd*', '*gd*'}) ;

% define list of structures
structure_list = {'alveus', 'ca_gm', 'ca_wm', 'subiculum', 'sub_wm', 'dg'}  ;

% subsample
for structure = structure_list
    %-- subsample volume
    aniso_seg.(char(structure)).mat = in_seg.(char(structure)).mat(:, 1:subsample_factor:end, :) ;
    %-- update matrix dimensions
    aniso_seg.(char(structure)).dim_mat = size(aniso_seg.(char(structure)).mat) ;
    %-- update voxel size
    aniso_seg.(char(structure)).vox_size = in_seg.(char(structure)).vox_size ;
    aniso_seg.(char(structure)).vox_size(2) = aniso_seg.(char(structure)).vox_size(2) * subsample_factor ;
    %-- retain data type
    aniso_seg.(char(structure)).data_type = in_seg.(char(structure)).data_type ;
end

% save all subsampled segmentations in output directory
for structure = structure_list
    aniso_structure_filename = strcat(output_segmentations_path, '/', structure{1}, '.ima') ;
    save_ima(aniso_structure_filename, aniso_seg.(char(structure)), 'S16','b') ;
end
end
