function compute_avgthickness_volume(input_segmentations_path, input_map_path, joined_subject_list, result_prefix, output_computation_path, streamlengths_path)
% compute_avgthickness_volume For each hemisphere (left, right) of each
% subject (controls and TLE patients), compute:
% 1. volume of the hippocampal ribbon
% 2. average thickness at the central surface extracted from the
% hippocampal ribbon
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
%   - input_map_path (str): path to the folder that contains the central
%         surfaces along with associated thickness maps of all subjects
%         Used to compute average thicknesses of central surfaces / maps
%         extracted with the RKHS method.
%   - joined_subject_list (str): joined ('[...] [...] ... [...]') list
%         of the subjects in the group study
%   - output_computation_path (str): path to the .json output file where
%         the average thicknesses and ribbon will be stored
%   - result_prefix (str): prefix of the files output in the output
%         directory of maps computation.
%   - streamlengths_path (str): (optional) path to a folder containing
%         all streamlengths maps (for all subjects/sides). Used to
%         compute average thicknesses of central surfaces / maps
%         extracted with the Laplacian method.

% check optional arguments
streamlengths_provided = false ;
if nargin >= 6 
    streamlengths_provided = true ;
end

% read list of subjects
subject_list = strsplit(joined_subject_list, ' ') ;

% go through the list of subject
for subject = subject_list
    % go through all hippocampal sides
    for side = {'left', 'right'}
        % check subject type (control or patient)
        subject_type = get_subject_type(char(subject)) ;
        % compute ribbon volume
        %-- define path to the input segmentations
        subjectside_segmentation_path = strcat(input_segmentations_path, '/', subject_type, 's', '/', char(subject), '/', char(side)) ;
        %-- read input segmentation
        [ca_gm] = loadHippoStructure(subjectside_segmentation_path, {'*CA_GM*', '*ca_gm*', '*CA_SP*', '*ca_sp*'}) ;
        [subiculum] = loadHippoStructure(subjectside_segmentation_path, {'*SUBICULUM*', '*Subiculum*', '*subiculum*'}) ;
        %-- compute ribbon
        ribbon_mat = ca_gm.mat + subiculum.mat ;
        %-- get volume
        %---- ribbon voxel number
        ribbon_voxel_number = sum(ribbon_mat(:)) ;
        %---- multiply by voxel dimensions
        voxel_volume = ca_gm.vox_size(1)*ca_gm.vox_size(2)*ca_gm.vox_size(3) ;
        ribbon_volume = ribbon_voxel_number*voxel_volume ;
        %-- add to list of volumes
        output_computation.('volume').(char(subject)).(char(side)) = ribbon_volume ;
        % compute ribbon average thickness
        %-- check which type of thickness maps was provided as an input
        %   - if computation for RKHS maps, should be a path to folder
        %       containing subfolders [subject_type]/[subject_id]/[side]
        %       with file [result_prefix]_thick_[subject_id]_[side].ima
        %   - if computation for Laplacian maps, should be a path to a
        %       folder containing all files
        %       [result_prefix]_[subject_id]_[side]_streamlengths_1.mat
        %       Each of this file contains T1 and T2, which are overall
        %       lengths for the 'ascending' and 'descending' streamlines
        %       respectively
        %-- define path to the data produced while computing the
        % thickness map 
        if streamlengths_provided == false
            % no streamlenghts provided: RKHS maps. Use [...]thick[...].ima
            subjectside_thickmap_data_path = strcat(input_map_path, '/', subject_type, 's', '/', char(subject), '/', char(side), '/', result_prefix, '_', char(subject), '_', char(side), '.mat') ;
            subjectside_thickmap_data = load(subjectside_thickmap_data_path) ;
            %-- get minimum voxel size (used to scale the thickness values
            % stored in the thickness map image)
            min_vox_size = min(subjectside_thickmap_data.metaData.vox_size) ;
            %-- define path to the thickness map
            subjectside_thickmap_path = strcat(input_map_path, '/', subject_type, 's', '/', char(subject), '/', char(side), '/', result_prefix, '_', char(subject), '_', char(side), '_thick', '.ima') ;
            %-- read thickness map image
            subjectside_thickmap_vol = load_ima(subjectside_thickmap_path) ;
            %-- scale
            subjectside_thickmap_vol.mat = min_vox_size * subjectside_thickmap_vol.mat ;
            %-- get average thickness
            ribbon_avgthick = mean(subjectside_thickmap_vol.mat(subjectside_thickmap_vol.mat~=0)) ;
        else
            % streamlengths provided: Laplacian maps. Use [...]_streamlenghts_1.mat
            %-- get path to the streamlengths .mat file for current subject/side
            subjectside_streamlengths_path = strcat(streamlengths_path, '/', result_prefix, '_', subject, '_', side, '_streamlengths_1.mat') ;
            %-- read streamlengths
            %---- overall structure
            subjectside_streamlengths = load(subjectside_streamlengths_path{1}) ;
            %---- 'ascending' streamline (point to outer boundary)
            subjectside_T1 = subjectside_streamlengths.T1 ;
            %---- 'descending' streamline (point to inner boundary)
            subjectside_T2 = subjectside_streamlengths.T2 ;
            %-- deduce thickness (T1+T2)
            subjectside_thick_map = subjectside_T1 + subjectside_T2 ;
            %-- force NaN values to 0
            subjectside_thick_map(isnan(subjectside_thick_map)) = 0 ;
            %-- get average thickness
            ribbon_avgthick = mean(subjectside_thick_map(subjectside_thick_map~=0)) ;
        end
        %-- add to list of average thicknesses
        output_computation.('avgthickness').(char(subject)).(char(side)) = ribbon_avgthick ;

    end
end

% store volumes / average thicknesseses for each subject/side in .json file
disp(' ')
disp('store volume / avg. thickness to .json file') ;
write_text_to_file(jsonencode(output_computation), output_computation_path) ;
disp(['volume / avg. thickness stored in', output_computation_path])
end


function subject_type = get_subject_type(subject)
% get_subject_type(subject). Get subject type ('control', 'patient')
% from subject ID
%
% Arguments:
%    - subject (str): subject ID (e.g., 'control1', 'patient4', etc.)
%
% Returns:
%    - subject_type (str): 'control' or 'patient'

% check if subject is a control subject

if ~isempty(strfind(subject, 'control'))
    subject_type='control' ;
elseif ~isempty(strfind(subject, 'patient'))
    subject_type='patient' ;
else
    disp('')
    disp('')
    disp('')
    disp('Error: Subject ID should contain control or patient.') ;
    disp('')
    disp('')
    disp('')
    exit
end
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
