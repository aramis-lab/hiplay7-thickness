#! /usr/bin/python

"""Code for Experiment III. B.

Script to run the experiments reported in section III. 'Experiments and
Results', subsection B. 'Influence of inter-rater variability' of
article 'A Diffeomorphic Vector Field Approaach to Analyze the Thickness
of the Hippocampus from 7T MRI' (Guyot et al.)
This runs in turn:
- part 4) 'Central surface and thickness estimation': computation on the
left and right hippocampus of each of the four control subjects of a
central surface with associated thickness map for a kernel size=10.
- part 5) 'Inter-rater analysis': 1. Inter-rater quantitative comparison
of each segment of each of the eight hippocampi (left/right) of the four
control subjects between two raters. 2. Average of the statistics across
the eight hippocampi.
- part 6) 'Comparison to Laplacian': Compute central surfaces /
thickness maps using the Laplacian method and perform the inter-rater
analysis using the Laplacian-derived surfaces.
"""

import os
import errno
import subprocess
import json
import argparse

import scipy as sp
import scipy.io
import joblib

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
    cli_description = 'Launcher for experiment B.'
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


def part4(
        code_path, input_data_path, output_data_path,
        subject_list, side_list, rater_list, n_cores):
    """part 4) 'Central surface and thickness estimation'

    Computes the central surface and thickness map obtained for a kernel
    size of 10 on all hippocampi for Experiment B.: left and right
    hippocampi of subjects 1, 2, 3 and 4.

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection B.
        output_data_path (str): path to the output data for 
            experiment III., subsection B.
        subject_list (list of str): list of the four subjects
        side_list (list of str): list of the two hippocampus sides
            'left' and 'right'
        rater_list (list of str): list of the two raters
        n_cores (int): >0. number of cores used for the parallel
            execution overs subject/side/rater

    Returns:
        N/A
    """
    # define input/output directories
    #-- input
    segmentations_input_data_path = os.path.join(
        input_data_path,
        "segmentations")
    #-- output
    surface_output_data_path = os.path.join(
        output_data_path,
        "4-centralSurfaceAndThicknessEstimation")
    os.makedirs(surface_output_data_path, exist_ok=True)

    # set kernel size to 10
    sigma = 10

    # run central surface / thickness code for sigma=10
    # done for each subject/side/rater
    in_tuple_list = []
    for subject in subject_list:
        for side in side_list:
            for rater in rater_list:
                #-- define the input segmentation path
                input_segmentation_path = os.path.join(
                    segmentations_input_data_path,
                    rater,
                    subject,
                    side)
                #-- create the output thickness dir 
                output_thickness_path = os.path.join(
                    surface_output_data_path,
                    rater,
                    subject,
                    side)
                os.makedirs(output_thickness_path, exist_ok=True)
                #-- add to input list
                result_prefix = 'hippo_thicknessMap'
                in_tuple = (
                    code_path,
                    input_segmentation_path,
                    output_thickness_path,
                    result_prefix,
                    sigma)
                in_tuple_list.append(in_tuple)
    #-- run function for all possible combinations of subject/side/rater
    joblib.Parallel(n_jobs=n_cores)(
        joblib.delayed(thick.compute_centralsurface_thicknessmap)(*in_tuple)
        for in_tuple in in_tuple_list)

    # save Matlab figures
    colourbar_min = 0.2
    colourbar_max = 2.3
    subject = subject_list[0]
    result_prefix = 'hippo_thicknessMap'
    in_tuple_list = []
    for side in side_list:
        for rater in rater_list:
            # define path to thickness map
            subside_output_data_path = os.path.join(
                surface_output_data_path,
                rater,
                subject,
                side)
            thickmap_path = os.path.join(
                subside_output_data_path,
                '{0}.mat'.format(result_prefix))
            # add to input list
            in_tuple = (
                code_path, thickmap_path, subside_output_data_path,
                result_prefix, colourbar_min, colourbar_max)
            in_tuple_list.append(in_tuple)
    #-- run function
    joblib.Parallel(n_jobs=n_cores)(
        joblib.delayed(thick.save_thickness_figure)(*in_tuple)
        for in_tuple in in_tuple_list)


def part5_compareindividual(
        code_path,
        input_data_path,
        output_data_path,
        subject_list,
        side_list,
        rater_list,
        result_prefix,
        n_cores):
    """part 5) 'Influence of kernel size'. Individual comparisons

    Compare quantitatively the central surface / thickness map obtained
    at sigma=10 to the other surfaces / maps obtained between rater 1
    and rater 2  using the three following metrics:
    1. correlation between thickness estimates
    2. mean distances between corresponding points
    3. mean absolute thickness differences.
    Surfaces obtained at different kernel sizes are matched by computing
    for each vertex of the source surface its closest vertex on the
    target surface.

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection B.
        output_data_path (str): path to the output data for 
            experiment III., subsection B.
        subject_list (list of str): list of the four subjects
        side_list (list of str): list of the two hippocampus sides
            'left' and 'right'
        rater_list (list of str): list of the two raters
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.
        n_cores (int): >0. number of cores used for the parallel
            execution overs subject/side/rater

    Returns:
        N/A
    """
    # define input/output directories
    #-- input
    boundaries_input_data_path = os.path.join(
        input_data_path,
        "boundaries")
    computemaps_output_data_path = os.path.join(
        output_data_path,
        "4-centralSurfaceAndThicknessEstimation")
    #-- output
    compareind_output_data_path = os.path.join(
        output_data_path,
        "5-interRaterAnalysis",
        "1-compareIndividual")
    os.makedirs(compareind_output_data_path, exist_ok=True)

    # Carry out surface quantitative comparison
    #-- define first and second raters
    rater1 = rater_list[0]
    rater2 = rater_list[1]
    #-- define list of all arguments for the above function
    in_tuple_list = []
    for subject in subject_list:
        for side in side_list:
            #-- check head/body/tail/whole boundaries in input folder
            boundaries_dict = None
            segment_boundaries_path = os.path.join(
                boundaries_input_data_path,
                "segment_boundaries_{0}_{1}.json".format(subject, side))
            with open(segment_boundaries_path, 'r') as boundaries_file:
                boundaries_dict = json.load(boundaries_file)
            if boundaries_dict is None:
                raise ValueError('No boundaries could be read')
            #-- create output compare directory
            individual_output_data_path = os.path.join(
                compareind_output_data_path,
                subject,
                side)
            os.makedirs(individual_output_data_path, exist_ok=True)
            #---- define reference / compared surfaces
            #------ reference surface: rater 1
            ref_surface_path = os.path.join(
                computemaps_output_data_path,
                rater1,
                subject,
                side,
                '{0}.mat'.format(result_prefix))
            #------ compared surface: rater 2
            compared_surface_path = os.path.join(
                computemaps_output_data_path,
                rater2,
                subject,
                side,
                '{0}.mat'.format(result_prefix))
            #------ output .json
            comparison_path = os.path.join(
                individual_output_data_path,
                'compare_individual.json')
            #------ output .mat (reference vertices and thicknesses
            # + corresponding compared vertices and thicknesses)
            comparison_data_path = os.path.join(
                individual_output_data_path,
                'compare_individual_data.mat')
            #-- define list of tuples
            noneIdentifier = None
            in_tuple = (
                ref_surface_path,
                compared_surface_path,
                segment_boundaries_path,
                comparison_path,
                comparison_data_path,
                noneIdentifier,
                code_path)
            in_tuple_list.append(in_tuple)
    #-- run function compute comparison metrics
    joblib.Parallel(n_jobs=n_cores)(
        joblib.delayed(thick.maps_quantitative_compare)(*in_tuple)
        for in_tuple in in_tuple_list)


def part5_compareaverage(
        input_data_path, output_data_path, subject_list, side_list):
    """part 5) 'Inter-rater analysis'. Average comparisons

    Average the quantitative comparisons computed on each left/right
    hippocampi of the four subjects (where the comparisons were carried
    out between the central surface / thickness maps obtained from the
    segmentations manually done by rater 1 and rater 2).

    Args:
        input_data_path (str): path to the input data for 
            experiment III., subsection B.
        output_data_path (str): path to the output data for 
            experiment III., subsection B.
        subject_list (list of str): list of the four subjects
        side_list (list of str): list of the two hippocampus sides
            'left' and 'right'

    Returns:
        N/A
    """
    # define input/output directories
    #-- input
    comparison_data_path = os.path.join(
        output_data_path,
        "5-interRaterAnalysis",
        "1-compareIndividual")
    #-- output
    compareavg_output_data_path = os.path.join(
        output_data_path,
        "5-interRaterAnalysis",
        "2-compareAverage")
    os.makedirs(compareavg_output_data_path, exist_ok=True)

    # initialise sum comparison 'left'/'right' dictionaries
    average_comp_dict = dict()
    sum_comp_dict = dict()
    metric_list = [
        'thickness_correlation',
        'mean_abs_thickness_difference',
        'mean_inter_surface_distance']
    segment_list = ['head', 'body', 'tail', 'whole']
    for side in side_list:
        sum_comp_dict[side] = dict()
        for metric in metric_list:
            sum_comp_dict[side][metric] = dict()
            for segment in segment_list:
                sum_comp_dict[side][metric][segment] = 0

    # initialise hippocampus count to 0
    hippocampus_count_dict = dict()
    for side in side_list:
        hippocampus_count_dict[side] = 0

    # read all the individual comparisons
    for subject in subject_list:
        for side in side_list:
            #-- define path to individual comparison
            individual_comp_path = os.path.join(
                comparison_data_path,
                subject,
                side,
                'compare_individual.json')
            #-- read individual comparison
            with open(individual_comp_path, 'r') as individual_comp_file:
                individual_comp_dict = json.load(individual_comp_file)
            #-- append to the sum of comparison metrics
            for metric in metric_list:
                for segment in segment_list:
                    sum_comp_dict[side][metric][segment] += individual_comp_dict[metric][segment]
            #-- increase hippocampus_count
            hippocampus_count_dict[side] += 1

    # average all the metrics by dividing by the hippocampus count
    average_comp_dict = dict()
    for side in side_list:
        average_comp_dict[side] = dict()
        for metric in metric_list:
            average_comp_dict[side][metric] = dict()
            for segment in segment_list:
                average_comp_dict[side][metric][segment] = 1.0*sum_comp_dict[side][metric][segment]/hippocampus_count_dict[side]

    # store the average comparison
    for side in side_list:
        average_comparison_path = os.path.join(
            compareavg_output_data_path,
            'compare_average_{0}.json'.format(side))
        with open(average_comparison_path, 'w') as average_comparison_file:
            json.dump(average_comp_dict[side], average_comparison_file)


def part5(
        code_path, input_data_path, output_data_path,
        subject_list, side_list, rater_list, n_cores):
    """part 5) 'Inter-rater analysis'

    Perform quantitative comparison of the central surface / thickness
    map obtained between the first rater and the second rater using the
    three following metrics:
    1. correlation between thickness estimates
    2. mean distances between corresponding points
    3. mean absolute thickness differences.
    The comparisons are first done for all subjects / hippocampus sides.
    They are then averaged to get a final output on each hippocampal
    segment.
    Surfaces obtained at different kernel sizes are matched by computing
    for each vertex of the source surface its closest vertex on the
    target surface.

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection B.
        output_data_path (str): path to the output data for 
            experiment III., subsection B.
        subject_list (list of str): list of the four subjects
        side_list (list of str): list of the two hippocampus sides
            'left' and 'right'
        rater_list (list of str): list of the two raters
        n_cores (int): >0. number of cores used for the parallel
            execution overs subject/side/rater

    Returns:
        N/A
    """
    # define data common to computation of maps and their comparisons
    result_prefix = 'hippo_thicknessMap'

    # compute central surfaces and thickness maps for each individual
    # between rater 1 and rater 2
    part5_compareindividual(
        code_path,
        input_data_path,
        output_data_path,
        subject_list,
        side_list,
        rater_list,
        result_prefix,
        n_cores)

    # average comparisons across all subjects
    part5_compareaverage(
        input_data_path,
        output_data_path,
        subject_list,
        side_list,)


def main():
    """Main function.
    
    Run parts
    4) 'Central surface and thickness estimation' and
    5) 'Inter-rater analysis'
    of subsection B. 'Influence of inter-rater variability' of
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
    #-- define experiment B specific input folder
    input_data_path = os.path.join(
        base_input_data_path,
        "2_influenceOfInterRaterVariability",
        "inputData")
    #-- number of cores for parallel execution
    n_cores = 1
    if args.n_cores:
        if args.n_cores > 0:
            n_cores = args.n_cores

    # handle I/O
    new_input_data_path = thick.process_io(input_data_path, output_data_path)

    # define list of subjects, hippocampus sides, raters
    subject_list = ['subject1', 'subject2', 'subject3', 'subject4']
    side_list = ['left', 'right']
    rater_list = ['rater1', 'rater2']

    # run part 4) 'Central surface and thickness estimation'
    part4(
        code_path, new_input_data_path, output_data_path,
        subject_list,side_list, rater_list, n_cores)

    # run part 5) 'Inter-rater analysis'
    part5(
        code_path, new_input_data_path, output_data_path,
        subject_list, side_list, rater_list, n_cores)


if __name__ == "__main__":
    main()
