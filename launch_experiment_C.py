#! /usr/bin/python

"""Code for Experiment III. C.

Script to run the experiments reported in section III. 'Experiments and
Results', subsection C. 'Application to in vivo 7T MRI group studies' of
article 'A Diffeomorphic Vector Field Approaach to Analyze the Thickness
of the Hippocampus from 7T MRI' (Guyot et al.)
This runs:
- part 4) 'Group study': 1. computation of the central surfaces and
associated thickness maps for each  patient and  control with a kernel
size = 10, 2. Computation of template for group 1 'controls +
contralateral patients' and group 2 'controls + ipsilateral patients',
3. Projection of the thickness maps of all controls and patients onto
their corresponding template, 4. Computation of average thickness for
each projected thickness map, 5. Spearman's rank correlation coefficient
comparison between average thickness and the volume of the hippocampal
ribbon, 6. Computation of central surfaces / thickness maps using the
Laplace method, 7. Computation of average thickness for the Laplace
derived surfaces, 8. Computation of Spearman's rank correlation
coefficients between average, 9. Computation of effect sizes between
controls and contralateral patients, as well as effect sizes between
controls and ipsilateral patients for volumes, RKHS thicknesses and
Laplace thicknesses on both left and right hemispheres.
"""

import os
import shutil
import errno
import subprocess
import json
import argparse
import collections

import numpy as np
import scipy.stats
import scipy.io
import joblib

import centralSurfacesAndThicknessCode.python_utils.matlab_run as mlabrun
import centralSurfacesAndThicknessCode.python_utils.thickness_functions as thick


def read_cli_args():
    """Read command-line interface arguments

    Parse the input to the command line with the argparse module.

    Args:
        N/A

    Returns:
        args (argparse.Namespace): parsed arguments
    """
    # read command line arguments
    cli_description = 'Launcher for experiment C.'
    parser = argparse.ArgumentParser(description=cli_description)
    # add arguments
    #-- mandatory arguments
    #---- input data
    parser.add_argument(
        'input_folder',
        help='path to input folder')
    #---- output data
    parser.add_argument(
        'output_folder',
        help='path to output folder')
    #-- optional arguments
    #---- number of cores
    parser.add_argument(
        '-n',
        '--n_cores',
        help='number of cores for parallel execution',
        type=int)
    # parse all arguments
    args = parser.parse_args()

    return args


def get_subside_result_prefix(result_prefix, subject_name, side):
    """Return a result prefix based on subject ID and side

    Used to force output central surfaces to have a unique name, so they
    can later all be moved to the same folder containing a data .xml
    parameter file to be input to the Deformetrica software.

    Args:
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.
        subject_name (str): subject ID
        side (str): hippocampus side. 'left' or 'right'
    """
    subside_result_prefix = str.format(
            '{0}_{1}_{2}', result_prefix, subject_name, side)

    return subside_result_prefix


def part4_centralsurfaceandthicknessestimation(
        code_path,
        input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list,
        result_prefix,
        n_cores):
    """part 4) 'Group study'. Computation of central surface / thickness

    Computes the central surface and thickness map obtained for a kernel
    size of 10 on all hippocampi for Experiment C.: left and right
    hippocampi of all controls and all patients

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection C.
        output_data_path (str): path to the output data for 
            experiment III., subsection C.
        control_list (list of str): list of the nine control subjects
        patient_list_dict2 (dictionary of dictionaries): list of the
            patients depending on epilepsy type and hippocampus side
            (patient_list_dict2[hippocampus_side][group])
        side_list (list of str): list of the two hippocampus sides
            'left' and 'right'
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.
        n_cores (int): >0. number of cores used for the parallel
            execution overs subject/side

    Returns:
        N/A
    """
    # define input/output directories
    #-- output
    surface_output_data_path = os.path.join(
        output_data_path,
        "4-groupStudy",
        "1-centralSurfaceAndThicknessEstimation")
    os.makedirs(surface_output_data_path, exist_ok=True)

    # assemble list of patients
    patient_list = (
        patient_list_dict2['left']['ipsi']
        + patient_list_dict2['left']['contra'])

    # define list of subjects based on the list of controls and the list
    # of patients. Each element defined as (subject_name, subject_type).
    Subject = collections.namedtuple('Subject', 'name type')
    subject_list = (
        [Subject(control, 'control') for control in control_list]
        + [Subject(patient, 'patient') for patient in patient_list])

    # set kernel size to 10
    sigma = 10

    # run central surface / thickness code for sigma=10
    # done for each subject/side/rater
    in_tuple_list = []
    for subject in subject_list:
        for side in side_list:
            #-- define the input segmentation path
            input_segmentation_path = os.path.join(
                input_data_path,
                '{0}s'.format(subject.type),
                subject.name,
                side)
            #-- create the output thickness dir 
            output_thickness_path = os.path.join(
                surface_output_data_path,
                '{0}s'.format(subject.type),
                subject.name,
                side)
            os.makedirs(output_thickness_path, exist_ok=True)
            #-- run the central surface / thickness code
            subside_result_prefix = get_subside_result_prefix(
                result_prefix, subject.name, side)
            in_tuple = (
                code_path,
                input_segmentation_path,
                output_thickness_path,
                subside_result_prefix,
                sigma)
            in_tuple_list.append(in_tuple)
    #-- run function for all possible combinations of subjects/side
    joblib.Parallel(n_jobs=n_cores)(
        joblib.delayed(thick.compute_centralsurface_thicknessmap)(*in_tuple)
        for in_tuple in in_tuple_list)


def sidegroup_to_objid(side, group):
    """Get object ID from both side and group.

    Convenience function to be used in both the generation of template
    for groups of controls and patients in left and right hemisphere and
    in the projection of surfaces onto the template.

    Args:
        side (str): 'left' or 'right'
        group (str): 'ipsi' (controls + ipsilateral TLE patients) or
            'contra' (controls + contralateral TLE patients)

    Returns:
        obj_id (str): object ID from side and group
    """
    objid = "{0}_{1}".format(side, group)

    return objid


def part4_centralsurfacetemplates(
        code_path,
        input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list,
        group_list,
        result_prefix):
    """part 4) 'Group study'. Computation of central surface templates

    Computes the central surface templates on left and right sides for
    two groups:
    - group 1: 'controls + ipsilateral patients'
    - group 2: 'controls + contralateral patients'

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection C.
        output_data_path (str): path to the output data for 
            experiment III., subsection C.
        control_list (list of str): list of the nine control subjects
        patient_list_dict2 (dictionary of dictionaries): list of the
            patients depending on epilepsy type and hippocampus side
            (patient_list_dict2[hippocampus_side][group])
        side_list (list of str): list of the two hippocampus sides
            'left' and 'right'
        group_list (list of str): list of epilepy type ('ipsilateral',
            and 'contralateral')
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.

    Returns:
        N/A
    """
    # define input/output directories
    #-- input
    computemaps_data_path = os.path.join(
        output_data_path,
        "4-groupStudy",
        "1-centralSurfaceAndThicknessEstimation")
    #-- intermediary
    intermediary_data_path = os.path.join(
        output_data_path,
        "debug",
        "4-groupStudy",
        "2-centralSurfaceTemplates")
    os.makedirs(intermediary_data_path, exist_ok=True)
    #-- output
    template_output_data_path = os.path.join(
        output_data_path,
        "4-groupStudy",
        "2-centralSurfaceTemplates")
    os.makedirs(template_output_data_path, exist_ok=True)

    # create input to the template creation software for each side/group
    # will generate three parameter files: optimisation, model and data
    #-- create intermediary folders
    for side in side_list:
        for group in group_list:
            os.makedirs(
                os.path.join(intermediary_data_path, side, group),
                exist_ok=True)
    #-- optimisation parameter : common to all subjects
    optimisation_path = os.path.join(
        intermediary_data_path,
        "optimisation_parameters.xml")
    thick.build_optimisation_file(optimisation_path)
    #-- model
    #---- model parameters
    template_sigma = 1.0
    data_kernel_width = 3
    deformation_kernel_width = 5
    #---- loop through all sides / groups
    objid_dict2 = dict()
    model_path_dict2 = dict()
    for side in side_list:
        model_path_dict2[side] = dict()
        objid_dict2[side] = dict()
        for group in group_list:
            #------ from first central surface, compute input
            # initialisation to the central surface software
            side_group_patient0 = patient_list_dict2[side][group][0]
            subside_result_prefix = get_subside_result_prefix(
                result_prefix, side_group_patient0, side)
            first_central_surface_path = os.path.join(
                computemaps_data_path,
                'patients',
                side_group_patient0,
                side,
                '{0}.vtk'.format(subside_result_prefix))
            init_template_path = os.path.join(
                intermediary_data_path,
                side,
                group,
                'init_template.vtk')
            thick.create_initial_template(
                first_central_surface_path,
                init_template_path)
            #------ define object ID corresponding to side and group
            objid_dict2[side][group] = sidegroup_to_objid(side, group)
            #------ generate the .xml file
            model_path_dict2[side][group] = os.path.join(
                intermediary_data_path,
                side,
                group,
                "model_parameters.xml")
            thick.build_model_file(
                model_path_dict2[side][group],
                'DeterministicAtlas',
                objid_dict2[side][group],
                init_template_path,
                template_sigma,
                data_kernel_width,
                deformation_kernel_width)
    #-- data
    data_path_dict2 = dict()
    subject_surfaces_dict2 = dict()
    dataparam_subject_surfaces_dict2 = dict()
    for side in side_list:
        data_path_dict2[side] = dict()
        subject_surfaces_dict2[side] = dict()
        dataparam_subject_surfaces_dict2[side] = dict()
        #---- define list of central surfaces for each controls for the
        # current side
        control_surface_list = [
            os.path.join(
                computemaps_data_path,
                'controls',
                control,
                side,
                str.format(
                    '{0}.vtk',
                    get_subside_result_prefix(result_prefix, control, side))
                ) for control in control_list]
        #---- get list of central surfaces for each patient group and merge
        # with controls
        for group in group_list:
            group_surface_list = [
                os.path.join(
                    computemaps_data_path,
                    'patients',
                    patient,
                    side,
                    str.format(
                        '{0}.vtk',
                        get_subside_result_prefix(result_prefix, patient, side))
                    ) for patient in patient_list_dict2[side][group]]
            subject_surfaces_dict2[side][group] = (
                control_surface_list
                + group_surface_list)
            #---- move all the input surfaces to the same directory as
            # data_parameters.xml (as required by Deformetrica)
            data_path_dict2[side][group] = os.path.join(
                intermediary_data_path,
                side,
                group,
                'data_parameters.xml')
            dataparam_subject_surfaces_dict2[side][group] = list()
            for subject_surface in subject_surfaces_dict2[side][group]:
                dataparam_subject_surface = os.path.join(
                    os.path.dirname(data_path_dict2[side][group]),
                    os.path.basename(subject_surface))
                shutil.copy(
                    subject_surface,
                    dataparam_subject_surface)
                dataparam_subject_surfaces_dict2[side][group].append(
                    dataparam_subject_surface)
        #---- write data file
        for group in group_list:
            subject_list = control_list + patient_list_dict2[side][group]
            thick.build_data_file(
                data_path_dict2[side][group],
                dataparam_subject_surfaces_dict2[side][group],
                subject_list,
                objid_dict2[side][group])

    # run the template creation software (Deformetrica) for each side/group
    for side in side_list:
        for group in group_list:
            #-- define output folder
            template_output_folder_path = os.path.join(
                template_output_data_path,
                side,
                group)
            template_output_log_path = os.path.join(
                template_output_folder_path,
                'deformetrica.log')
            template_output_errlog_path = os.path.join(
                template_output_folder_path,
                'deformetrica.err.log')
            #-- create folder
            os.makedirs(template_output_folder_path, exist_ok=True)
            #-- build command
            template_create_cmd = str.format(
                'deformetrica estimate {0} {1} -p {2} --output={3} > {4} 2> {5}',
                model_path_dict2[side][group],
                data_path_dict2[side][group],
                optimisation_path,
                template_output_folder_path,
                template_output_log_path,
                template_output_errlog_path)
            #-- run command
            print('Run')
            print(template_create_cmd)
            subprocess.run(template_create_cmd, shell=True)


def part4_templateprojections(
        code_path,
        input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list,
        group_list,
        result_prefix):
    """part 4) 'Group study'. Projection onto central surface templates

    For each side/group (e.g., 'left'/'ipsi'), select all associated
    surfaces and project them onto the central surface template
    previously obtained with Deformetrica.

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection C.
        output_data_path (str): path to the output data for 
            experiment III., subsection C.
        control_list (list of str): list of the nine control subjects
        patient_list_dict2 (dictionary of dictionaries): list of the
            patients depending on epilepsy type and hippocampus side
            (patient_list_dict2[hippocampus_side][group])
        side_list (list of str): list of the two hippocampus sides
            'left' and 'right'
        group_list (list of str): list of epilepy type ('ipsilateral',
            and 'contralateral')
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.

    Returns:
        N/A
    """
    # define input/output directories
    #-- input
    computemaps_data_path = os.path.join(
        output_data_path,
        "4-groupStudy",
        "1-centralSurfaceAndThicknessEstimation")
    computetemplates_data_path = os.path.join(
        output_data_path,
        "4-groupStudy",
        "2-centralSurfaceTemplates")
    #-- intermediary
    intermediary_data_path = os.path.join(
        output_data_path,
        "debug",
        "4-groupStudy",
        "2-centralSurfaceTemplates")
    os.makedirs(intermediary_data_path, exist_ok=True)
    #-- output
    projection_output_data_path = os.path.join(
        output_data_path,
        "4-groupStudy",
        "3-templateProjections")
    os.makedirs(projection_output_data_path, exist_ok=True)

    # create output folders for each group in the group study
    for side in side_list:
        for group in group_list:
            #-- define output folder
            sidegroup_output_data_path = os.path.join(
                projection_output_data_path,
                side,
                group)
            #-- create folder
            os.makedirs(sidegroup_output_data_path, exist_ok=True)

    # for each group in the study, project thickness onto template
    for side in side_list:
        for group in group_list:
            # define 
            # Create median surface within Matlab
            #-- define list of paths required by Matlab
            matlab_path_list = [
                os.path.join(code_path, "core"),
                os.path.join(code_path, "externalLibraries", "matlab")]
            #--  define .m function to be called
            matlab_func = 'computeAvgThickness_views'
            #-- define list of all arguments for the above function
            #---- final template
            final_template_filename = str.format(
                'DeterministicAtlas__EstimatedParameters__Template_{0}.vtk',
                sidegroup_to_objid(side, group))
            final_template_path = os.path.join(
                computetemplates_data_path,
                side,
                group,
                final_template_filename)
            #---- initial template: dummy, not used
            initial_template_path = 'dummy'
            #---- output data folder
            sidegroup_output_data_path = os.path.join(
                projection_output_data_path,
                side,
                group)
            #---- list of subjects for the current group
            subject_list = control_list + patient_list_dict2[side][group]
            #---- list of initially computed subject surfaces
            #------ controls
            init_control_surf_path_list = [
                os.path.join(
                    computemaps_data_path,
                    'controls',
                    control,
                    side,
                    '{0}_{1}_{2}.vtk'.format(result_prefix, control, side)
                ) for control in control_list]
            #------ patients
            init_patient_surf_path_list = [
                os.path.join(
                    computemaps_data_path,
                    'patients',
                    patient,
                    side,
                    '{0}_{1}_{2}.vtk'.format(result_prefix, patient, side)
                    ) for patient in patient_list_dict2[side][group]]
            #------ assemble controls and patients
            init_sub_surf_path_list = (
                init_control_surf_path_list
                + init_patient_surf_path_list)
            join_init_sub_surf_path_list = str.format(
                '\'{0}\'',
                ' '.join(init_sub_surf_path_list))
            #---- list of corresponding reconstructed subject surfaces
            sidegroup_computetemplates_data_path = os.path.join(
                computetemplates_data_path,
                side,
                group)
            templ_sub_surf_filename_list = [
                str.format(
                    'DeterministicAtlas__Reconstruction__{0}__subject_{1}.vtk',
                    sidegroup_to_objid(side, group),
                    subject) for subject in subject_list]
            templ_sub_surf_path_list = [
                os.path.join(
                    sidegroup_computetemplates_data_path,
                    templ_sub_surf_filename
                ) for templ_sub_surf_filename in templ_sub_surf_filename_list]
            join_templ_sub_surf_path_list = str.format(
                '\'{0}\'',
                ' '.join(templ_sub_surf_path_list))
            #---- output thickness file prefix
            output_file_prefix = 'thickness_projection'
            #---- test title: dummy, not used
            test_title = 'dummy'

            matlab_arg_list = [
                final_template_path,
                initial_template_path,
                sidegroup_output_data_path,
                join_init_sub_surf_path_list,
                join_templ_sub_surf_path_list,
                output_file_prefix,
                test_title]
            #-- run in Matlab
            mlabrun.run_matlab_func(
                matlab_func, matlab_arg_list, matlab_path_list)


def part4_avgThicknessVolumeComputation(
        code_path,
        input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list,
        result_prefix):
    """part 4) 'Group study'. Computation of avg. thickness and volume

    For each subject (both controls and TLE ipsi/contra -lateral
    patients) in the group study, each side, compute:
    - the volume of the hippocampal ribbon (CA-SP + subiculum)
    - the average thickness found on the hippocampal ribbon central
        surface extracted at a previous stage

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection C.
        output_data_path (str): path to the output data for 
            experiment III., subsection C.
        control_list (list of str): list of the nine control subjects
        patient_list_dict2 (dictionary of dictionaries): list of the
            patients depending on epilepsy type and hippocampus side
            (patient_list_dict2[hippocampus_side][group])
        side_list (list of str): list of the two hippocampus sides
            'left' and 'right'
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.

    Returns:
        N/A
    """
    # define input/output directories
    #-- input
    segmentation_data_path = os.path.join(
        input_data_path)
    computemaps_data_path = os.path.join(
        output_data_path,
        "4-groupStudy",
        "1-centralSurfaceAndThicknessEstimation")
    #-- output
    avgthickvol_output_data_path = os.path.join(
        output_data_path,
        "4-groupStudy",
        "4-avgThicknessVolumeComputation")
    os.makedirs(avgthickvol_output_data_path, exist_ok=True)

    # assemble list of patients
    patient_list = (
        patient_list_dict2['left']['ipsi']
        + patient_list_dict2['left']['contra'])

    # define list of subjects based on the list of controls and the list
    # of patients. Each element defined as (subject_name, subject_type).
    Subject = collections.namedtuple('Subject', 'name type')
    subject_list = control_list + patient_list
    joined_subject_list = str.format(
        '\'{0}\'',
        ' '.join(subject_list))

    # define path to .json to store volumes and average thicknesses
    avgthickness_volume_path = os.path.join(
        avgthickvol_output_data_path,
        'avgthickness_volume.json')

    # run Matlab code
    #-- define function name
    matlab_func = 'compute_avgthickness_volume'
    #-- define list of arguments
    matlab_arg_list = [
        segmentation_data_path,
        computemaps_data_path,
        joined_subject_list,
        result_prefix,
        avgthickness_volume_path]
    #-- define paths to be added to the Matlab interpreter
    matlab_path_list = [
        os.path.join(code_path, "core"),
        os.path.join(code_path, "externalLibraries", "matlab"),
        os.path.join(code_path, "extra") ]
    #-- run in Matlab
    mlabrun.run_matlab_func(
        matlab_func, matlab_arg_list, matlab_path_list)

    # plot graphs with Matplotlib
    thick.plot_avgthick_volume(
        avgthickness_volume_path,
        control_list,
        patient_list_dict2,
        side_list,
        avgthickvol_output_data_path,
        volume_label = '(b) Volume (mm$^3$)',
        avgthick_label = '(a) Average Thickness (mm)',
        avgthick_legend=False)


def part4_spearmanCorrelationsComputation(
        input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list):
    """part 4) 'Group study'. Computation of Spearman correlations

    Compute Spearman correlation between the hippocampal ribbon volume
    and average thickness of the corresponding central surface for each
    subject of the group study in both left and right hemispheres.

    Args:
        input_data_path (str): path to the input data for 
            experiment III., subsection C.
        output_data_path (str): path to the output data for 
            experiment III., subsection C.
        control_list (list of str): list of the nine control subjects
        patient_list_dict2 (dictionary of dictionaries): list of the
            patients depending on epilepsy type and hippocampus side
            (patient_list_dict2[hippocampus_side][group])
        side_list (list of str): list of the two hippocampus sides
            'left' and 'right'
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.

    Returns:
        N/A
    """
    # define input/output directories
    #-- input
    computemeasures_data_path = os.path.join(
        output_data_path,
        "4-groupStudy",
        "4-avgThicknessVolumeComputation")
    #-- output
    spearman_output_data_path = os.path.join(
        output_data_path,
        "4-groupStudy",
        "5-spearmanCorrelationsComputation")
    os.makedirs(spearman_output_data_path, exist_ok=True)

    # retrieve volume and average thickness for each patient/side
    avgthickness_volume_path = os.path.join(
        computemeasures_data_path,
        "avgthickness_volume.json")
    avgthick_volume_dict = None
    with open(avgthickness_volume_path, 'r') as avgthickness_volume_file:
        avgthick_volume_dict = json.load(avgthickness_volume_file)
    if avgthick_volume_dict is None:
        raise ValueError(
            'No hippocampal ribbon volume / average thickness found')

    # initialise Spearman correlation dictionary (left and right side)
    spearman_corr_dict = dict()

    # go through all sides
    for side in side_list:
        # get list of volumes / average thicknesses for all subjects
        #-- initialise lists
        volume_list = []
        avgthick_list = []
        #-- go through all subjects
        for subject in avgthick_volume_dict['volume']:
            volume_list.append(
                avgthick_volume_dict['volume'][subject][side])
            avgthick_list.append(
                avgthick_volume_dict['avgthickness'][subject][side])
        # compute Spearman correlation
        spearman_corr_dict[side] = scipy.stats.spearmanr(
                avgthick_list, volume_list)[0]

    # store Spearman correlation dictionary in output .json file
    spearman_corr_path = os.path.join(
        spearman_output_data_path,
        'spearman_correlations.json')
    with open(spearman_corr_path, 'w') as spearman_corr_file:
        json.dump(spearman_corr_dict, spearman_corr_file)


def part4(
        code_path,
        input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list,
        group_list,
        n_cores):
    """part 4) 'Group study'

    Conduct all the experiments reported in part 4 'Group study' of
    section C. 'Application to in vivo 7T MRI group studies':
    1. computation of the central surfaces and associated thickness maps
    for each  patient and  control with a kernel size = 10
    2. Computation of template for group 1 'controls + ipsilateral
    patients' and group 2 'controls + contralateral patients'
    3. Projection of the thickness maps of all controls and patients
    onto their corresponding template
    4. Computation of average thickness for each projected thickness map
    5. Spearman's rank correlation coefficient comparison between
    average thickness and the volume of the hippocampal ribbon.

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection C.
        output_data_path (str): path to the output data for 
            experiment III., subsection C.
        control_list (list of str): list of the nine control subjects
        patient_list_dict2 (dictionary of dictionaries): list of the
            patients depending on epilepsy type and hippocampus side
            (patient_list_dict2[hippocampus_side][group])
        side_list (list of str): list of the two hippocampus sides
            'left' and 'right'
        group_list (list of str): list of epilepy type ('ipsilateral',
            and 'contralateral')
        n_cores (int): >0. number of cores used for the parallel
            execution overs subject/side

    Returns:
        N/A
    """
    # define data common to computation of maps and their comparisons
    result_prefix = 'hippo_thicknessMap'

    # Compute central surfaces and associated thickness maps for all
    # controls and patients
    part4_centralsurfaceandthicknessestimation(
        code_path,
        input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list,
        result_prefix,
        n_cores)

    # compute templates central surface for
    # - group 1: 'controls + ipsilateral patients'
    # - group 2: 'controls + contralateral patients'
    part4_centralsurfacetemplates(
        code_path,
        input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list,
        group_list,
        result_prefix)

    # For each side/group, project each individual thickness
    # map onto the corresponding side/group template
    part4_templateprojections(
        code_path,
        input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list,
        group_list,
        result_prefix)

    # For each side, each subject, compute average thickness and volume
    part4_avgThicknessVolumeComputation(
        code_path,
        input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list,
        result_prefix)

    # For each side, each subject, compute average thickness and volume
    part4_spearmanCorrelationsComputation(
        input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list)


def main():
    """Main function.
    
    Run parts
    4) 'Group study',
    of subsection C. 'Application to in vivo 7T MRI group studies' of
    section III. 'Experiments and Results'

    Args:
        N/A
Returns:
        N/A
    """
    # retrieve path to code
    launcher_path = os.path.abspath(__file__)
    experiments_path = os.path.dirname(launcher_path)
    code_path = os.path.join(
        experiments_path,
        "centralSurfacesAndThicknessCode")

    # read command-line arguments
    #-- paths to base input / output folder
    args = read_cli_args()
    base_input_data_path = args.input_folder
    output_data_path = args.output_folder
    #-- define experiment C specific input folder
    input_data_path = os.path.join(
        base_input_data_path,
        "3_applicationToInVivo7tMriGroupStudy",
        "inputData")
    #-- number of cores for parallel execution
    n_cores = 1
    if args.n_cores:
        if args.n_cores > 0:
            n_cores = args.n_cores

    # handle I/O
    new_input_data_path = thick.process_io(input_data_path, output_data_path)

    # define lists
    #-- controls
    control_list = [
        'control1',
        'control2',
        'control3',
        'control4',
        'control6',
        'control8',
        'control9',
        'control10',
        'control11']
    #-- patients groups
    group_list = ['ipsi', 'contra']
    patient_list_dict2 = dict()
    #---- left side
    patient_list_dict2['left'] = dict()
    patient_list_dict2['left']['ipsi'] = [
        'patient3',
        'patient5',
        'patient8',
        'patient9',
        'patient12']
    patient_list_dict2['left']['contra'] = [
        'patient1',
        'patient2',
        'patient6']
    #---- right side
    patient_list_dict2['right'] = dict()
    patient_list_dict2['right']['ipsi'] = patient_list_dict2['left']['contra']
    patient_list_dict2['right']['contra'] = patient_list_dict2['left']['ipsi']
    #-- hippocampus sides
    side_list = ['left', 'right']

    # run part 4) 'Group study'
    part4(
        code_path,
        new_input_data_path,
        output_data_path,
        control_list,
        patient_list_dict2,
        side_list,
        group_list,
        n_cores)


if __name__ == "__main__":
    main()
