function maps_quantitative_compare(ref_map_path, comp_map_path, segment_boundaries_path, identifier, comparison_path, comparison_data_path, oversample)
% maps_quantitative_compare Quantitatively compare two
% {central surfaces + thickness maps} derived from segmentations of a
% ribbon of hippocampal gray matter
%
% Arguments:
%   - ref_map_path (str): path to the
%         {central surfaces + thickness maps} structure that contains the
%         reference surface
%   - comp_map_path (str): path to the
%         {central surfaces + thickness maps} structure that contains the
%         surface to be compared to the reference surface
%   - segment_boundaries_path (str): path to the .json file that stores a 
%         structure containing the start and end positions of head, body
%         tail and whole for both reference and compared maps.
%         The boundaries have been visually selected from the
%         RIA-oriented volume using Freeview
%   - identifier (str): identifies which segment is being compared
%   - comparison_path (str): path to the .json file where all the results
%         of the quantitative comparison will be stored
%   - comparison_data_path (str): path to the .mat file to store
%         intermediary data used to generate the final values for each
%         of the comparison metric
%   - oversample (int): identifies if oversample option is provided
%         legacy, not in use any longer

% percentile set to 100, do not threshold any thing
prctile_value = 100

% check optional arguments
oversample_factor = 1 ;
if nargin >= 7
    % oversample factor
    % use when the images have been oversampled to make them isotropic
    % (Laplacian experiment for the 0.2 x 0.3 x 0.2 post-mortem specimen)
    oversample_factor = str2num(oversample) ;
end

% load input data
disp(' ')
disp('load data: surfaces, segment boundaries')
%-- {central surfaces + thickness maps} from reference and compared maps
ref_map = load(ref_map_path) ;
comp_map = load(comp_map_path) ;
%-- segment boundaries
segment_boundaries = jsondecode(fileread(segment_boundaries_path)) ;


% for each point of the compared central surface, compute its closest
% counterpart on the reference surface, as well as the distance between
% the two points
% Note: For this to be accurate, all points on the compared central
%     surface must have a counterpart on the reference surface
%     (e.g., if we are downsampling an isotropic surface to simulate
%     anisotropy, holes will likely appear in the anisotropic surface,
%     meaning that we cannot use the isotropic surface as a reference,
%     as some of its points correspond to the holes in the simulated
%     anisotropic image. Instead, we have to choose the anisotropic
%     surface as the reference and the isotropic surface as the
%     compared surface)
disp(' ')
disp('compute vertex correspondences between surfaces')
[ref_map_closest,ref_map_dist] = dsearchn(ref_map.SkelVert,comp_map.SkelVert);


% retrieve thickness for each closest point on the compared surface
disp(' ')
disp('compute corresponding thicknesses')
ref_map_thick = ref_map.SkelThick(ref_map_closest);

% remove NaN values from thickness values
comp_map_thick = comp_map.SkelThick ;


% deduce difference in thickness between all points of the reference
% surface and their counterpart on the compared surface
diff_thick = (ref_map_thick-comp_map_thick);

% save differences from subsampled to non-sample to file
%-- get path of folders containing all the maps
[map_folder, dummy, dummy] = fileparts(ref_map_path) ;
%-- define new diff structure containing central surface information
diffmap.SkelVert = comp_map.SkelVert ;
diffmap.SkelFace = comp_map.SkelFace ;
diffmap.SkelThick = -diff_thick ;
%-- save the figure
if strcmp(identifier, 'interrater')
    % interrater: colour in
    colourbar_min = -1 ;
    colourbar_max = 2.3 ;
else
    colourbar_min = -1 ;
    colourbar_max = 2.3 ;
end
thickfig_prefix = strcat('diff_', identifier) ;
save_thicknessmap_figure(diffmap, map_folder, thickfig_prefix, colourbar_min, colourbar_max) ;



% subdivide the central surface into segments head, body, tail, whole:
%-- skel vertices ([map].Skel) have their positions in the same space
% as the cropped volume within which the central surface and thickness
% values were extracted. However, the start/end boundaries for each
% segment were retrieved within 1. with Freeview, 2. within the whole
% volume. We therefore will
% 1. reverse the boundaries to reflect the fact that Brainvisa reverts
%    the Y coordinates compared to those read in Freeview
% 2. subtract the origin of the cropped volume from the boundaries
% remove the gap between the origin (0,0,0) of the whole volume and
% the origin of the cropped volume to the start/end boundaries, so
% that they match the space of {the cropped volume} / {[map].Skel}
disp(' ')
disp('split central surfaces into segments')
whole_n = double(ref_map.metaData.dim_mat(2)) ;
origin_n = double(ref_map.metaData.origin(2)) ;
%-- create binary masks for each segments (i.e., vertices of the
% central surface that belongs to the segment are labelled 1, while
% all the other vertices are labelled 0)
% The boundaries are computed in the 'mm' coordinate system from the
% reference image. For The computation to be accurate, the slice
% boundaries given in the input segment_boundaries_path file must have
% been retrieved from the reference image
for segment = {'head', 'body', 'tail', 'whole'}
    %---- initialise by labelling all vertices as belonging to segment
    subdiv_map.(char(segment)) = ones(size(comp_map.SkelVert,1),1) ;
    %---- get boundaries for current segment
    %------ voxel boundaries in cropped volume
    vox_crop_seg_start = (whole_n -segment_boundaries.(char(segment)).end) - origin_n ;
    vox_crop_seg_end = (whole_n - segment_boundaries.(char(segment)).start) - origin_n ;
    %------ mm boundaries in cropped volume
    % (see calcule_skel_3D(...) to undertand how mm to vox is computed)
    min_vox_size = min(ref_map.metaData.vox_size) ;
    aniso_y = min_vox_size/ref_map.metaData.vox_size(2) ;
    mm_crop_seg_start = (vox_crop_seg_start-1.0)/aniso_y ;
    mm_crop_seg_end = (vox_crop_seg_end-1.0)/aniso_y ;
    %------ mm boundaries in overall volume (common to both ref and
    % compared surfaces)
    % (see hippocamp_Seg2Skel(...) to understand how cropped mm to
    % overall mm is computed)
    ref_origin = double(ref_map.metaData.origin) ;
    ref_vox_size = double(ref_map.metaData.vox_size) ;
    mm_seg_start = (ref_origin(2)-1)*ref_vox_size(2)+mm_crop_seg_start*min_vox_size*oversample_factor ;
    mm_seg_end = (ref_origin(2)-1)*ref_vox_size(2)+mm_crop_seg_end*min_vox_size*oversample_factor ;
    %---- exclude vertices (x,y,z) where y<start and y>end
    % Note: the subdivision is conducted for the compared image, as we
    %     wish to to get the point correspondences from compared to
    %     reference and not vice versa (see previous comment)
    %------ y<start
    subdiv_map.(char(segment))(comp_map.SkelVert(:,2)<mm_seg_start) = 0 ;
    %------ y>end
    subdiv_map.(char(segment))(comp_map.SkelVert(:,2)>mm_seg_end) = 0 ;
end

% add elements to the comparison data
%-- vertices, distances and thicknesses
comparison_data.ref_map_SkelVert = ref_map.SkelVert ; 
comparison_data.comp_map_SkelVert = comp_map.SkelVert ; 
comparison_data.ref_map_closest = ref_map_closest ; 
comparison_data.ref_map_dist = ref_map_dist ; 
comparison_data.ref_map_thick = ref_map_thick ; 
comparison_data.comp_map_thick = comp_map_thick ; 
comparison_data.diff_thick = diff_thick ; 
%-- head/body/tail subdivision for all vertices of compared map
comparison_data.subdiv_map = subdiv_map ; 

% carry out quantitative comparisons of surface for all the segments
disp(' ')
disp('quantitatively compare reference surface to the other surface:') ;
%-- define the names of the different comparison methods, as will appear
% in the file storing all the results
corr_method = 'thickness_correlation' ;
thickdiff = 'mean_abs_thickness_difference' ;
surfdiff = 'mean_inter_surface_distance' ;
%-- compute the comparison measures and store into structure
for segment = {'head', 'body', 'tail', 'whole'}
    %---- retrieve mask corresponding to the current segment
    segment_mask = logical(subdiv_map.(char(segment))) ;
    %---- compute mean absolute thickness differences between the two
    % surfaces
    %------ remove nan values
    nan_segment_diff_thick = diff_thick(segment_mask) ;
    segment_diff_thick = nan_segment_diff_thick(~isnan(nan_segment_diff_thick)) ;
    %------ only retain vertices for which the absolute thickness
    % difference falls below the percentile threshold
    prctile_thresh = prctile(abs(segment_diff_thick), prctile_value) ;
    prctile_thick_vertices = (abs(segment_diff_thick) <= prctile_thresh) ;
    segment_diff_thick = segment_diff_thick(prctile_thick_vertices) ;
    %---- compute without the nan values
    comparison.(char(thickdiff)).(char(segment)) = mean(abs(segment_diff_thick)) ;
    %---- compute correlation between reference and compared surfaces
    %------ remove NaN values
    nan_segment_ref_map_thick = ref_map_thick(segment_mask) ;
    nan_segment_comp_map_thick = comp_map_thick(segment_mask) ;
    segment_ref_map_thick = nan_segment_ref_map_thick(~isnan(nan_segment_ref_map_thick.*nan_segment_comp_map_thick)) ;
    segment_comp_map_thick = nan_segment_comp_map_thick(~isnan(nan_segment_ref_map_thick.*nan_segment_comp_map_thick)) ;
    %------ only retain vertices for which the absolute thickness
    % difference falls below the percentile threshold
    segment_ref_map_thick = segment_ref_map_thick(prctile_thick_vertices) ;
    segment_comp_map_thick = segment_comp_map_thick(prctile_thick_vertices) ;
    %------ compute without the nan values
    comparison.(char(corr_method)).(char(segment)) = corr(segment_ref_map_thick, segment_comp_map_thick) ;
    %---- compute mean inter-surface distance between the two surfaces
    segment_ref_map_dist = ref_map_dist(segment_mask) ;
    comparison.(char(surfdiff)).(char(segment)) = mean(abs(segment_ref_map_dist)) ;
    %---- disp
    disp(' ')
    disp(['segment: ', char(segment)]) ;
    disp(['thickness correlation: ', num2str(comparison.(char(corr_method)).(char(segment)))]) ;
    disp(['mean abs. thickness difference: ', num2str(comparison.(char(thickdiff)).(char(segment)))]) ;
    disp(['mean inter-surface distance: ', num2str(comparison.(char(surfdiff)).(char(segment)))]) ;
    %---- add elements to the comparison data
    %------ segment mask
    comparison_data.segments.(char(segment)).segment_mask = segment_mask ;
    %------ thickness correlations
    comparison_data.segments.(char(segment)).segment_ref_map_thick = segment_ref_map_thick ;
    comparison_data.segments.(char(segment)).segment_comp_map_thick = segment_comp_map_thick ;
    comparison_data.segments.(char(segment)).(char(corr_method)) = comparison.(char(corr_method)).(char(segment)) ;
    %------ mean absolute thickness differences
    comparison_data.segments.(char(segment)).diff_thick = segment_diff_thick ;
    comparison_data.segments.(char(segment)).(char(thickdiff)) = comparison.(char(thickdiff)).(char(segment)) ;
    %------ mean inter surface distances
    comparison_data.segments.(char(segment)).segment_ref_map_dist = segment_ref_map_dist ; 
    comparison_data.segments.(char(segment)).(char(surfdiff)) = comparison.(char(surfdiff)).(char(segment)) ; 
end

% store all quantitative comparisons to .json file
disp(' ')
disp('store comparisons to .json file') ;
write_text_to_file(jsonencode(comparison), comparison_path) ;
disp(['comparisons stored in', comparison_path])

% store intermediary data to file
save(comparison_data_path, 'comparison_data') ;

end



function write_text_to_file(my_str, input_file_path)
% write_text_to_file Write a string to a user-defined file
%
% Arguments:
%   - my_str (str): string to be written in text file
%   - input_file_path (str): path to the file where the string
%       will be stored
fid = fopen(input_file_path, 'wt') ;
fprintf(fid, '%s\n', my_str) ;
fclose(fid) ;
end
