#! /usr/bin/python

""" Code for Experiment III. A.

Script to run the experiments reported in section III. 'Experiments and
Results', subsection A. 'Robustness with respect to kernel size and
anisotropy' of article 'A Diffeomorphic Vector Field Approaach to
Analyze the Thickness of the Hippocampus from 7T MRI' (Guyot et al.)
This runs in turn:
- part 4) 'Central surface and thickness estimation': computation of a
central surface with associated thickness map for a kernel size=10
- part 5) 'Influence of kernel size': computation of the central
surfaces and their associated thickness maps for kernel sizes in
[3, 5, 10, 15], then quantitative comparison of the central surfaces /
thickness maps obtained at a kernel size of 10 to the other
surfaces / maps obtained for all the other kernel sizes.
- part 6) 'Influence of anisotropy': generation of pseudo anisotropic
segmentations, extraction of central surfaces / thickness maps and
quantitative comparison to the original central surface / thickness map.
- part 7) 'Comparison to Laplacian': extracts central surfaces /
thickness maps for all subsampling factor and perform quantitative
comparison of surfaces at subsampling factor >1 to the central surface /
thickness map at subsampling factor =1.
"""

import os
import errno
import subprocess
import json
import argparse

import scipy as sp
import scipy.io

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
    cli_description = 'Launcher for experiment A.'
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
    # parse all arguments
    args = parser.parse_args()

    return args


def part4(code_path, input_data_path, output_data_path):
    """part 4) 'Central surface and thickness estimation'

    Computes the central surface and thickness map obtained for a kernel
    size of 10.

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection A.
        output_data_path (str): path to the output data for 
            experiment III., subsection A.

    Returns:
        N/A
    """
    # define input/output directories
    #-- input
    segmentation_input_data_path = os.path.join(
        input_data_path,
        "segmentation")
    #-- output
    surface_output_data_path = os.path.join(
        output_data_path,
        "4-centralSurfaceAndThicknessEstimation")
    os.makedirs(surface_output_data_path, exist_ok=True)

    # set kernel size to 10
    sigma = 10

    # run central surface / thickness code for sigma=10
    result_prefix = 'hippo_thicknessMap_{0}'.format(sigma)
    thick.compute_centralsurface_thicknessmap(
        code_path,
        segmentation_input_data_path,
        surface_output_data_path,
        result_prefix,
        sigma)
    
    # Save Matlab figures
    thickmap_path = os.path.join(
        surface_output_data_path,
        '{0}.mat'.format(result_prefix))
    colourbar_min = 0.2
    colourbar_max = 2.3
    thick.save_thickness_figure(
        code_path, thickmap_path, surface_output_data_path, result_prefix,
        colourbar_min, colourbar_max)


def part5_computemaps(
        code_path, input_data_path, output_data_path, result_prefix):
    """part 5) 'Influence of kernel size'. Compute maps.

    Computes the central surfaces and thickness maps for kernel sizes in
    [3, 5, 10, 15].

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection A.
        output_data_path (str): path to the output data for 
            experiment III., subsection A.
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.

    Returns:
        N/A
    """
    # define input/output directories
    #-- input
    segmentation_input_data_path = os.path.join(
        input_data_path,
        "segmentation")
    #-- output
    maps_output_data_path = os.path.join(
        output_data_path,
        "5-influenceOfKernelSize",
        "1-computeMaps")
    os.makedirs(maps_output_data_path, exist_ok=True)

    # define list of possible kernel sizes
    sigma_list = [3, 5, 10, 15]

    # run central surface / thickness code for all defined kernels
    for sigma in sigma_list:
        # append kernel size to result prefix
        sigma_result_prefix = '{0}_{1}'.format(result_prefix, sigma)
        thick.compute_centralsurface_thicknessmap(
            code_path,
            segmentation_input_data_path,
            maps_output_data_path,
            sigma_result_prefix,
            sigma)


def part5_comparemaps(
        code_path, input_data_path, output_data_path, result_prefix):
    """part 5) 'Influence of kernel size'

    Compare quantitatively the central surface / thickness map obtained
    at sigma=10 to the other surfaces / maps obtained for all the other
    kernel sizes [3, 5, 15] using the three following metrics:
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
            experiment III., subsection A.
        output_data_path (str): path to the output data for 
            experiment III., subsection A.
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.

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
        "5-influenceOfKernelSize",
        "1-computeMaps")
    #-- output
    comparemaps_output_data_path = os.path.join(
        output_data_path,
        "5-influenceOfKernelSize",
        "2-compareMaps")
    os.makedirs(comparemaps_output_data_path, exist_ok=True)

    # check head/body/tail/whole boundaries in input folder
    boundaries_dict = None
    segment_boundaries_path = os.path.join(
        boundaries_input_data_path, "segment_boundaries.json")
    with open(segment_boundaries_path, 'r') as boundaries_file:
        boundaries_dict = json.load(boundaries_file)
    if boundaries_dict is None:
        raise ValueError('No boundaries could be read')

    # Carry out surface quantitative comparison
    #-- define list of all arguments for the above function
    #---- reference surface: kernel size 10
    ref_surface_path = os.path.join(
        computemaps_output_data_path,
        result_prefix+'_10.mat')
    #---- compared surface: all kernel sized in [3, 5, 15]
    sigma_list = [3, 5, 15]
    for sigma in sigma_list:
        #------ compared surface
        compared_surface_path = os.path.join(
            computemaps_output_data_path,
            '{0}_{1}.mat'.format(result_prefix, sigma))
        #------ output .json
        comparison_path = os.path.join(
            comparemaps_output_data_path,
            'compare_sigma10_to_sigma{0}.json'.format(sigma))
        #------ output .mat (reference vertices and thicknesses
        # + corresponding compared vertices and thicknesses)
        comparison_data_path = os.path.join(
            comparemaps_output_data_path,
            'compare_original_to_sigma{0}_data.mat'.format(sigma))
        #-- compute comparison metrics
        thick.maps_quantitative_compare(
            ref_surface_path,
            compared_surface_path,
            segment_boundaries_path,
            comparison_path,
            comparison_data_path,
            sigma,
            code_path)



def part5(code_path, input_data_path, output_data_path):
    """part 5) 'Influence of kernel size'

    Computes the central surfaces and thickness maps for kernel sizes in
    [3, 5, 10, 15]. Compare quantitatively the 
    central surface / thickness map obtained at sigma=10 to the other
    surfaces / maps obtained for all the other kernel sizes using the
    three following metrics:
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
            experiment III., subsection A.
        output_data_path (str): path to the output data for 
            experiment III., subsection A.

    Returns:
        N/A
    """
    # define data common to computation of maps and their comparisons
    result_prefix = 'hippo_thicknessMap'

    # compute central surfaces and thickness maps for all kernels
    part5_computemaps(
        code_path, input_data_path, output_data_path, result_prefix)

    # compare surface;map between kernel size 10 and all other kernels
    part5_comparemaps(
        code_path, input_data_path, output_data_path, result_prefix)


def part6_createanisotropicsegmentations(
        code_path, input_data_path, output_data_path, result_prefix,
        subsample_list):
    """part 6) 'Influence of anisotropy'. Create aniso- segmentations

    Generate 'simulated' segmentations from anisotropic acquisitions by
    sub-sampling segmentations from a quasi-isotropic segmentation.

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection A.
        output_data_path (str): path to the output data for 
            experiment III., subsection A.
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.
        subsample_list (list): list of integers. List of all subsample
            factors, i.e., list of 'f' values, where each 'simulated'
            anisotropic segmentation is generated by taking one slice
            out of 'f' from the original 'quasi-isotropic' segmentation)

    Returns:
        N/A
    """
    # define input/output directories
    segmentation_input_data_path = os.path.join(
        input_data_path,
        "segmentation")
    aniso_output_data_path = os.path.join(
        output_data_path,
        "6-influenceOfAnisotropy",
        "1-createAnisotropicSegmentations")
    os.makedirs(aniso_output_data_path, exist_ok=True)

    # generate 'simulated' segmentation for all the subsample factors
    # (include subsampling=1 - no subsampling - in the computed maps)
    for subsample in [1]+subsample_list:
        #-- create the output segmentation dir 
        output_segmentation_path = os.path.join(
            aniso_output_data_path,
            'subsample_{0}'.format(subsample))
        os.makedirs(output_segmentation_path, exist_ok=True)
        #-- generate segmentation within output segmentation dir using
        # Matlab
        #---- define list of paths required by Matlab
        matlab_path_list = [
            os.path.join(code_path, "core"),
            os.path.join(code_path, "externalLibraries", "matlab"),
            os.path.join(code_path, "extra")]
        #----  define .m function to be called
        matlab_func = 'create_anisotropic_segmentations'
        #---- define list of all arguments for the above function
        matlab_arg_list = [
            segmentation_input_data_path,
            output_segmentation_path,
            subsample]
        #---- run in Matlab
        mlabrun.run_matlab_func(
            matlab_func, matlab_arg_list, matlab_path_list)


def get_subsample_list():
    """Return subsampling list for parts 6 and 7

    Each subsample value in the list corresponds to an anisotropic
    factor for an MR 'pseudo'-segmentation of the hippocampal ribbon.
    Thise defines a list of 'f' values, where each 'simulated'
    anisotropic segmentation is generated by taking one slice out of
    'f' from the original 'quasi-isotropic' segmentation.
    Using a function here, so this list of subsampling factors is
    accessible from both part 6 (robustness to anisotropy) and part 7
    (comparison to Laplacian potential method).
    """
    subsample_list = [2, 3, 4, 5, 6]

    return subsample_list


def part6_computemaps(
        code_path, input_data_path, output_data_path, result_prefix,
        subsample_list):
    """part 6) 'Influence of anisotropy'. Compute maps.

    Computes the central surfaces and thickness maps for all defined
    subsampling factors.

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection A.
        output_data_path (str): path to the output data for 
            experiment III., subsection A.
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.
        subsample_list (list): list of integers. List of all subsample
            factors, i.e., list of 'f' values, where each 'simulated'
            anisotropic segmentation is generated by taking one slice
            out of 'f' from the original 'quasi-isotropic' segmentation)

    Returns:
        N/A
    """
    # define output dir
    comparemaps_output_data_path = os.path.join(
        output_data_path,
        "6-influenceOfAnisotropy",
        "2-computeMaps")
    os.makedirs(comparemaps_output_data_path, exist_ok=True)

    # run central surface / thickness code for all subsample factors
    # (include subsampling=1 - no subsampling - in the computed maps)
    for subsample in [1]+subsample_list:
        #-- define input segmentation dir (dependent on subsampling)
        subsample_input_data_path = os.path.join(
            output_data_path,
            "6-influenceOfAnisotropy",
            "1-createAnisotropicSegmentations",
            "subsample_{0}".format(subsample))
        #-- append subsampling to result prefix
        subsample_result_prefix = '{0}_{1}'.format(result_prefix, subsample)
        #-- define kernel size sigma = 10
        sigma = 10
        #-- compute maps
        thick.compute_centralsurface_thicknessmap(
            code_path,
            subsample_input_data_path,
            comparemaps_output_data_path,
            subsample_result_prefix,
            sigma)

    # save Matlab figures
    colourbar_min = 0.2
    colourbar_max = 2.3
    for subsample in subsample_list:
        # define path to thickness map
        subsample_result_prefix = '{0}_{1}'.format(result_prefix, subsample)
        thickmap_path = os.path.join(
            comparemaps_output_data_path,
            '{0}.mat'.format(subsample_result_prefix)
        )
        thick.save_thickness_figure(
            code_path, thickmap_path, comparemaps_output_data_path,
            subsample_result_prefix,
            colourbar_min, colourbar_max)


def part6_comparemaps(
        code_path, input_data_path, output_data_path, result_prefix,
        subsample_list):
    """part 6) 'Influence of anisotropy'

    Compare quantitatively the central surface / thickness map obtained
    at for the reference segmentation to the other surfaces / maps
    obtained for all 'simulated' segmentations from anisotropic MR data
    using different subsample factors.
    The comparisons are achieved with the three following metrics:
    1. correlation between thickness estimates
    2. mean distances between corresponding points
    3. mean absolute thickness differences.
    Surfaces obtained for different anisotropy factor are matched by
    computing for each vertex of the source surface its closest vertex
    on the target surface.

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection A.
        output_data_path (str): path to the output data for 
            experiment III., subsection A.
        result_prefix (str): prefix of the files output in the
            output directory of maps computation.
        subsample_list (list): list of integers. List of all subsample
            factors, i.e., list of 'f' values, where each 'simulated'
            anisotropic segmentation is generated by taking one slice
            out of 'f' from the original 'quasi-isotropic' segmentation)

    Returns:
        N/A
    """
    # define input/output directories
    #-- input
    boundaries_input_data_path = os.path.join(
        input_data_path,
        "boundaries")
    createsegs_output_data_path = os.path.join(
        output_data_path,
        "6-influenceOfAnisotropy",
        "1-createAnisotropicSegmentations")
    computemaps_output_data_path = os.path.join(
        output_data_path,
        "6-influenceOfAnisotropy",
        "2-computeMaps")
    #-- output
    comparemaps_output_data_path = os.path.join(
        output_data_path,
        "6-influenceOfAnisotropy",
        "3-compareMaps")
    os.makedirs(comparemaps_output_data_path, exist_ok=True)

    # check head/body/tail/whole boundaries in input folder
    boundaries_dict = None
    segment_boundaries_path = os.path.join(
        boundaries_input_data_path, "segment_boundaries.json")
    with open(segment_boundaries_path, 'r') as boundaries_file:
        boundaries_dict = json.load(boundaries_file)
    if boundaries_dict is None:
        raise ValueError('No boundaries could be read')

    # Carry out surface quantitative comparison
    #-- define list of all arguments for the above function
    #---- reference surface: subsample factor 1 (no subsampling)
    ref_surface_path = os.path.join(
        computemaps_output_data_path,
        result_prefix+'_1.mat')
    #---- compared surface: all subsample factors
    for subsample in subsample_list:
        print('SUBSAMPLE: {0}'.format(subsample))
        #------ compared surface
        compared_surface_path = os.path.join(
            computemaps_output_data_path,
            '{0}_{1}.mat'.format(result_prefix, subsample))
        #------ output .json
        comparison_path = os.path.join(
            comparemaps_output_data_path,
            'compare_original_to_subsampling{0}.json'.format(subsample))
        #------ output .mat (reference vertices and thicknesses
        # + corresponding compared vertices and thicknesses)
        comparison_data_path = os.path.join(
            comparemaps_output_data_path,
            'compare_original_to_subsampling{0}_data.mat'.format(subsample))
        #-- compute comparison metrics
        thick.maps_quantitative_compare(
            ref_surface_path,
            compared_surface_path,
            segment_boundaries_path,
            comparison_path,
            comparison_data_path,
            subsample,
            code_path)

    # Compute anisotropy factor for each subsample factors and output
    # to .json file
    #-- compute the anisotropy factor for the original segmentation
    # (subsample=1)
    #---- read voxel dimensions from .mat stored CA-SP segmentation map
    casp_mat_path = os.path.join(
        createsegs_output_data_path,
        'subsample_1',
        'CA_SP.mat')
    casp = sp.io.loadmat(casp_mat_path)
    original_vox_size = casp['CA_SP'][0][0][2][0]
    original_aniso_factor = 1.0*original_vox_size[1]/original_vox_size[0]
    #-- compute anisotropy factor for each subsample factor as
    # {original anisotropy factor} * {subsample factor}
    anisotropy_factor = dict()
    for subsample in subsample_list:
        anisotropy_factor_value = original_aniso_factor*subsample
        subsample_factor = {'anisotropy_factor': anisotropy_factor_value}
        anisotropy_factor[
                'subsample_factor={0}'.format(subsample)] = subsample_factor
    #-- save to .json file
    anisotropy_factor_path = os.path.join(
        comparemaps_output_data_path,
        'anisotropy_factors.json')
    with open(anisotropy_factor_path, 'w') as anisotropy_factor_file:
        json.dump(anisotropy_factor, anisotropy_factor_file)
    print('')
    print('saved anisotropy factors in {0}'.format(anisotropy_factor_path))


def part6(code_path, input_data_path, output_data_path):
    """part 6) 'Influence of anisotropy'

    Generate 'simulated' segmentations from anisotropic acquisitions by
    sub-sampling segmentations from a quasi-isotropic segmentation.
    Compare quantitatively the central surface / thickness map obtained
    from the quasi-isotropic segmentation to the maps derived from all
    the 'simulated' segmentations. Quantitative comparison is achieved
    using the same measures as those in part5(...).

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the input data for 
            experiment III., subsection A.
        output_data_path (str): path to the output data for 
            experiment III., subsection A.

    Returns:
        N/A
    """
    # define data common to computation of maps and their comparisons
    result_prefix = 'hippo_thicknessMap'

    # define subsampling list
    # (i.e., define list of 'f' values, where each 'simulated'
    # anisotropic segmentation is generated by taking one slice out
    # of 'f' from the original 'quasi-isotropic' segmentation)
    subsample_list = get_subsample_list()

    # generate simulated segmentations from anisotropic acquisitions
    # with different subsample factors (each subsample factor
    # corresponds to an anisotropy factor)
    part6_createanisotropicsegmentations(
        code_path, input_data_path, output_data_path, result_prefix,
        subsample_list)

    # compute central surfaces and thickness maps for all the simulated
    # segmentations
    part6_computemaps(
        code_path, input_data_path, output_data_path, result_prefix,
        subsample_list)

    # compare surface;map between the reference map and maps computed on
    # all the simulated segmentations
    part6_comparemaps(
        code_path, input_data_path, output_data_path, result_prefix,
        subsample_list)


def main():
    """Main function.
    
    Run parts
    4) 'Central surface and thickness estimation',
    5) 'Influence of kernel size' and
    6) 'Influence of anisotropy'
    of subsection A. 'Robustness with respect to kernel size and
    anisotropy' of section III. 'Experiments and Results'

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
    #-- define experiment A specific input folder
    input_data_path = os.path.join(
        base_input_data_path,
        "1_robustnessWithRespectToKernelSizeAndAnisotropy",
        "inputData")

    # handle I/O
    new_input_data_path = thick.process_io(input_data_path, output_data_path)

    # run part 4) 'Central surface and thickness estimation'
    part4(code_path, new_input_data_path, output_data_path)

    # run part 5) 'Influence of kernel size'
    part5(code_path, new_input_data_path, output_data_path)

    # run part 6) 'Influence of anisotropy'
    part6(code_path, new_input_data_path, output_data_path)


if __name__ == "__main__":
    main()
