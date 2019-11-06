#! /usr/bin/python

"""Central surface / thickness utilities

Code to generate and manipulate central surfaces and their associated
thickness maps from segmentation data within user-provided folders.
Wraps Matlab code.
"""

import os
import shutil
import pickle
import json

import math
import decimal
import numpy as np
import scipy as sp
import scipy.io
import scipy.ndimage
import scipy.interpolate
import scipy.integrate
import matplotlib.pyplot as plt
import matplotlib.lines
import nibabel as nib
import nilearn as nil
import nilearn.image
import vtk

import centralSurfacesAndThicknessCode.python_utils.matlab_run as mlabrun

plt.switch_backend('agg')


def copy_input_files(input_data_path, output_data_path):
    """Copy all input files to new location

    Create 'new' input data folder, with the same tree structures for
    the subfolders.
    This is because the Matlab code used for this experiments write data
    to the input folder, something we wish to avoid.
    Tried using symlinks, does not seem to work

    Args
        input_data_path (str): path to the input data folder containing
            boundaries and segmentation subfolders
        output_data_path (str): output path provided by the user for
            this experiment. This will be used to store the 'new' input data.

    Returns:
        new_input_data_path (str): path to the 'new' input data folder
            where all symlinks will be stored
    """
    # create intermediary folder
    intermediary_data_path = os.path.join(
        output_data_path,
        "debug")
    os.makedirs(intermediary_data_path, exist_ok=True)
    # create new input data path
    new_input_data_path = os.path.join(
        intermediary_data_path,
        "inputData")
    shutil.copytree(input_data_path, new_input_data_path)

    return new_input_data_path


def process_io(input_data_path, output_data_path):
    """ Process I/O

    Handles the inputs/outputs data folder for experiment A:
    Check that the input data folder provided by the user is valid, that
    the output data they provided is empty, so we do not accidentaly
    overwrite any data and compute copy input data to debug location
    inside the output data folder (because some of the Matlab functions
    will write within the input data folder).

    Args:
        input_data_path (str): path to the input data folder containing
            boundaries and segmentation subfolders
        output_data_path (str): output path provided by the user for
            this experiment. This will be used to store the 'new' input data.

    Returns:
        new_input_data_path (str): path to the 'new' input data folder
            where all symlinks will be stored
    """
    # check input data folder is valid
    if not os.path.isdir(input_data_path):
        raise OSError('{0} is not a valid folder.'.format(input_data_path))

    # only accept empty data folders
    if os.path.isdir(output_data_path):
        if os.listdir(output_data_path):
            raise OSError('Please Provide an empty output folder')

    # Matlab code can potentially write data into the input data folder
    # so we are creating a 'new' input data folder inside our output
    # data folder, with the same subfolders and files
    new_input_data_path = copy_input_files(input_data_path, output_data_path)

    return new_input_data_path


def compute_centralsurface_thicknessmap(
        code_path,
        input_data_path,
        output_data_path,
        result_prefix,
        sigma):
    """Compute central surface / thickness map

    Compute the central surface and thickness map of a hippocampus for a
    specified kernel size. The segmentation of the hippocampal gray matter
    is retrieved from a user-provided directory.

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        input_data_path (str): path to the folder containing the input
            segmentation
        output_data_path (str): path to the output folder to store the
            central surfaces and thickness maps. Will be populated with
            (having [xxx]=[result_prefix]):
            - [xxx].mat: all the matrices / structures generated while
                running the central surface / thickness map algorithm.
                Contains the vector field, the central surface
                (skeleton) mesh vertices and faces, streamlines and
                metadata
            - [xxx].mesh: brainvisa mesh of the hippocampal central
                surface and associated thickness map
            - [xxx].vtk: VTK mesh of the hippocampal central
            - [xxx]_thick.dim/ima: brainvisa image of the thickness
                map across the hippocampal ribbon of gray matter
            - [xxx]_pot.dim/ima: brainvisa image of the potential
                map associated to the thickness map
            - [xxx]_skel.dim/ima: brainvisa binary mask image of the
                central surface (skeleton) of the hippocampal
                ribbon. Output for visualisation purposes.
        result_prefix (str): prefix of the files output in the
            output directory.
        sigma (int): value of the kernel size used in the algorithm

    Returns:
        N/A
    """
    # Create median surface within Matlab
    #-- define list of paths required by Matlab
    matlab_path_list = [
        os.path.join(code_path, "core"),
        os.path.join(code_path, "externalLibraries", "matlab")]
    #--  define .m function to be called
    matlab_func = 'create_median_surface'
    #-- define list of all arguments for the above function
    dimension = '3D' # all experiments run with the 3D version
    input_prefix = 'dummy' # legacy parameter, not in use any longer
    matlab_arg_list = [
        dimension,
        input_data_path,
        input_prefix,
        output_data_path,
        result_prefix,
        sigma]
    #-- run in Matlab
    mlabrun.run_matlab_func(matlab_func, matlab_arg_list, matlab_path_list)


def save_thickness_figure(
        code_path, thickmap_path, thickfig_folder, thickfig_prefix,
        colourbar_min, colourbar_max):
    """Save thickness map visualisation as a .png image

    Save the thickness map visualisation from two different views as a
    .png image.
    Wrapper to a Matlab function.

    Args:
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        thickmap_path (str): path to the .mat structure containing the
            central surface volume
        thickfig_folder (str): path to the folder where to output .png
            visualisation
        thickfig_prefix (str): prefix for the name of the files that
            will be created inside thickfig_folder
        colourbar_min (float): min value for the colourbar
        colourbar_max (float): max value for the colourbar

    Returns:
        N/A
    """
    # Create median surface within Matlab
    #-- define list of paths required by Matlab
    matlab_path_list = [
        os.path.join(code_path, "core"),
        os.path.join(code_path, "externalLibraries", "matlab"),
        os.path.join(code_path, "extra")]
    #--  define .m function to be called
    matlab_func = 'save_thickness_figure'
    #-- define list of all arguments for the above function
    matlab_arg_list = [
        thickmap_path,
        thickfig_folder,
        thickfig_prefix,
        colourbar_min,
        colourbar_max ]
    #-- run in Matlab
    mlabrun.run_matlab_func(matlab_func, matlab_arg_list, matlab_path_list)


def maps_quantitative_compare(
        ref_surface_path,
        compared_surface_path,
        segment_boundaries_path,
        comparison_path,
        comparison_data_path,
        identifier,
        code_path,
        oversample=None):
    """Quantitative comparison of two central surface / thickness

    Quantitatively compare two {central surfaces + thickness maps}
    derived from segmentations of a ribbon of hippocampal gray matter
    Three comparison metrics are computed:
    1. thickness correlations between the two thickness maps
    2. Mean absolute thickness differences between the two thickness
        maps
    3. Mean inter surface distances between the two central surfaces.
    Points correspondences between the two surfaces are obtained by
    taking all the vertices of the compared surface, and for each
    vertex, retrieve the closest vertex on the reference surface.
    Wraps a Matlab function.
    
    Args:
        ref_surface_path (str): path to the
            {central surfaces + thickness maps} Matlab structure .mat
            file that contains the reference surface
        comp_map_path (str): path to the
            {central surfaces + thickness maps} Matlab structure .mat
            file that contains the surface to be compared to the
            reference surface
        segment_boundaries_path (str): path to the .json file that
            stores a structure containing the start and end positions of
            head, body tail and whole for both reference and compared
            maps.
            The boundaries have been visually selected from the
            RIA-oriented volume using Freeview
        comparison_path (str): path to the .json file where all the
            results of the quantitative comparison will be stored
        comparison_data_path (str): path to the .mat file to store
            intermediary data used to generate the final values for each
            of the comparison metric
        identifier (str): identifies which map is compared
        code_path (str): path to the Matlab code to generate central
            surfaces and thickness maps
        oversample (int): identifies if oversample option is
            provided. Use when conducting experiments where images have
            been oversampled (e.g. Laplacian experiment)

    Returns:
        N/A
    """
    # derive the list of paths required by Matlab to run the surface
    # comparison script from the path to the code folder
    matlab_path_list = [
        os.path.join(code_path, "core"),
        os.path.join(code_path, "extra")]

    # define the .m function that will be called by Matlab
    matlab_func = 'maps_quantitative_compare'

    # define the list or arguments that will be passed to this Matlab
    # function
    matlab_arg_list = [
        ref_surface_path,
        compared_surface_path,
        segment_boundaries_path,
        identifier,
        comparison_path,
        comparison_data_path]
    if oversample is not None:
        matlab_arg_list.append(oversample)

    # run the function with the arguments using Matlab
    mlabrun.run_matlab_func(matlab_func, matlab_arg_list, matlab_path_list)


def read_polydata(vtk_path):
    """Read polydata from .vtk file

    Simple polydata reader.

    Args:
        vtk_path (String): name of the .vtk file

    Returns:
        polydata (vtk polydata): polydata from the file
    """
    pdReader = vtk.vtkPolyDataReader()
    pdReader.SetFileName(vtk_path)
    pdReader.Update()
    polydata = pdReader.GetOutput()

    return polydata


def write_polydata(polydata_path, polydata):
    """Write polydata to .vtk file

    Simple polydata writer

    Args:
        polydata_path: (String)
        polydata: (vtk polydata)

    Returns:
        N/A
    """
    pdWriter = vtk.vtkPolyDataWriter()
    pdWriter.SetInputData(polydata)
    pdWriter.SetFileName(polydata_path)
    pdWriter.Write()


def create_initial_template(central_surface_path, init_template_path):
    """Create initialisation to template creation

    Produces the initial template file (input to deformetrica) from one
    of the central surfaces (usually the first patient in the list of
    shapes).
    Repairs the central surface provided by the user (i.e. fill the
    holes) and makes it smoother.

    Args:
        central_surface_path (python string): one of the central
            surfaces computed from a segmentation of the hippocampus.
        init_template_path (python string)

    Returns:
        N/A
    """
    # produce two smooth central surfaces
    #-- 1) smooth with hard boundaries
    #---- read central surface
    central_surface_polydata = read_polydata(central_surface_path)
    #---- fill holes
    holes_filter = vtk.vtkFillHolesFilter()
    holes_filter.SetInputData(central_surface_polydata)
    holes_filter.SetHoleSize(5.0)
    holes_filter.Update()
    central_surface_repaired = holes_filter.GetOutput()
    #---- make smoother
    smooth_filter = vtk.vtkSmoothPolyDataFilter()
    smooth_filter.SetInputData(central_surface_repaired)
    smooth_filter.SetRelaxationFactor(1.0)
    smooth_filter.SetNumberOfIterations(15)
    smooth_filter.SetEdgeAngle(150.0)
    smooth_filter.Update()
    #---- output to intial template file
    smooth_surface = smooth_filter.GetOutput()
    write_polydata(init_template_path, smooth_surface)


def build_optimisation_file(optimisation_path):
    """Build optimisation parameters .xml file

    Generate optimisation .xml file, input to the central surface
    template creation software (Deformetrica).

    Args:
        optimisation_path (string): path to the .xml optimisation file

    Returns:
        N/A
    """
    with open(optimisation_path, 'w') as optim_file:
        optim_file.write('<?xml version="1.0"?>\n')
        optim_file.write('\n')
        optim_file.write('<optimization-parameters>\n')
        optim_file.write('<optimization-method-type>GradientAscent'\
                         '</optimization-method-type>\n')
        optim_file.write('<initial-step-size>0.001</initial-step-size>\n')
        optim_file.write('<max-iterations>100</max-iterations>\n')
        optim_file.write('<number-of-processes>17</number-of-processes>\n')
        optim_file.write('<freeze-control-points>Off</freeze-control-points>\n')
        optim_file.write('<freeze-template>Off</freeze-template>\n')
        optim_file.write('</optimization-parameters>')


def build_model_file(model_path, model_type, object_id, input_template_filename, template_sigma, data_kernel_width, deformation_kernel_width):
    """Build model parameters .xml file

    Generate model .xml file., input to the central surface template
    creation software (Deformetrica).

    Args:
        model_path (string): path to the .xml model file
        model_type (string): deformation type. 'DeterministicAtlas' or
            'Registration'
        object_id (string): template identifier, later retrieved in .xml
            data file
        input_template_filename (string): initialisation for
            Deformetrica template creation software
        template_sigma (positive float): Lagrange multiplier
        data_kernel_width (float): weight to measure the fit between
            registered surface and target
        deformation_kernel_width (float): weight to measure the
            smoothness of the deformation

    Returns:
        N/A
    """
    model_file = open(model_path, 'w')
    with open(model_path, 'w') as model_file:
        model_file.write('<?xml version="1.0"?>\n')
        model_file.write(' \n')
        model_file.write('<model>\n')
        model_file.write('<model-type>{0}</model-type>\n'.format(model_type))
        model_file.write('<dimension>3</dimension>\n')
        model_file.write(' \n')
        model_file.write('<template>\n')
        model_file.write('    <object id="{0}">\n'.format(object_id))
        model_file.write('        <filename>{0}</filename>\n'.format(
            input_template_filename))
        model_file.write('        <deformable-object-type>'\
                         'SurfaceMesh</deformable-object-type>\n')
        model_file.write('        <attachment-type>'\
                         'varifold</attachment-type>\n')
        model_file.write('        <kernel-width>{0}</kernel-width>\n'.format(
            data_kernel_width))
        model_file.write('        <kernel-type>torch</kernel-type>\n')
        model_file.write('        <noise-std>{0}</noise-std>\n'.format(
            template_sigma))
        model_file.write('    </object>\n')
        model_file.write('</template>\n')
        model_file.write(' \n')
        model_file.write('<deformation-parameters>\n')
        model_file.write('    <kernel-width>{0}</kernel-width>\n'.format(
            deformation_kernel_width))
        model_file.write('    <kernel-type>torch</kernel-type>\n')
        model_file.write('    <number-of-timepoints>10'\
                         '</number-of-timepoints>\n')
        model_file.write('</deformation-parameters>\n')
        model_file.write('</model>')


def build_data_file(data_path, surface_list, subject_list, object_id):
    """Build data .xml data file

    Generate data .xml file, input to the template creation software
    (Deformetrica).

    Args:
        data_path (string)
        surface_list (list of string): list of all the surfaces required
            to build the template. 
        subject_list (list of string): list of all the corresponding
            subjects (controls and ipsi/contra -lateral TLE patients)
        object_id (string): identifier of the template, previously
            defined in the .xml model file

    Returns:
        N/A
    """
    # count number of surfaces
    subject_number = len(subject_list)
    # produce .xml from all the surfaces provided to the function.
    # Note: we are only using the base file name of each surface (as
    #       opposed to the entire path), because Deformetrica needs the
    #       surfaces to be in the same folder as the data file and file
    #       paths are generated from both the path to data.xml and each
    #       individual surface base file name. Moving the surfaces to
    #       the correct data.xml folder is done prior to running the
    #       current function
    with open(data_path, 'w') as data_file:
        data_file.write('<?xml version="1.0"?>\n')
        data_file.write('<data-set>\n')
        data_file.write('\n')
        for subject_index in range(subject_number):
            surface = surface_list[subject_index]
            subject = subject_list[subject_index]
            data_file.write('    <subject id="{0}">\n'.format(subject))
            data_file.write('        <visit id="experiment">\n')
            data_file.write('            <filename object_id="{0}">{1}</filename>\n'.format(
                object_id, os.path.basename(surface)))
            data_file.write('        </visit>\n')
            data_file.write('    </subject>\n')
        data_file.write('\n')
        data_file.write('\n')
        data_file.write('</data-set>')


def plot_avgthick_volume(
        avgthickness_volume_path,
        control_list,
        patient_list_dict2,
        side_list,
        output_data_path,
        volume_label=None,
        avgthick_label=None,
        avgthick_legend=False,
        laplacian=False):
    """Plot avg. thickness and volumes for all populations

     Will plot 1. the average thickness and 2. the volumes for controls,
     contralateral patients and ipsilateral patients for both left and
     right hemispheres.

     Args:
         avgthickness_volume_path (str): path to the .json file that
             contains the average thickness and volumes for each subject
             (control, contralateral patient or ipsilateral patient) for
             both left and right hemispheres
         control_list (list of str): list of the nine control subjects
         patient_list_dict2 (dictionary of dictionaries): list of the
             patients depending on epilepsy type and hippocampus side
             (patient_list_dict2[hippocampus_side][group])
         side_list (list of str): list of the two hippocampus sides
             'left' and 'right'
         output_data_path (str): folder where the two output graphs 
             (1. average thickness for both left and right hemispheres,
             2. volume for both left and right hemispheres)
         volume_label (str): label under the X axis for the volume plot
         avgthick_label (str): label under the X axis for the avg.
             thickness plot
         avgthick_legend (boolean): if True, add a legend to the avg.
            thickness plot, else do not add the legend.
            RKHS maps do not need this option, however it is useful for
            the avg. thicknesses obtained with Laplace's method.
        laplacian (boolean): indicates if we are  plotting volumes and
            average thicknesses for the RKHS method or for the Laplace
            method. Used to change the Y scale depending on which method
            was used.

     Returns:
         N/A
    """
    # read files containing all out measures
    with open(avgthickness_volume_path, 'r') as measure_dict:
        measure_dict = json.load(measure_dict)

    # read volumes and average thicknesss
    volume_dict =  measure_dict['volume']
    avgthick_dict =  measure_dict['avgthickness']

    # define colours
    colour_dict2 = dict()
    colour_dict2['left'] = dict()
    colour_dict2['left']['control'] = '#cbcb00'
    colour_dict2['left']['contra'] = '#03b700'
    colour_dict2['left']['ipsi'] = '#0093ce'
    colour_dict2['right'] = dict()
    colour_dict2['right']['control'] = '#8800ff'
    colour_dict2['right']['contra'] = '#ff00d3'
    colour_dict2['right']['ipsi'] = '#ff7f00'

    # generate plots for all measures and sides
    # in each case, plot controls, ipsi-lateral patients and
    # contra-lateral patients
    measure_list = ['volume', 'avgthickness']
    for measure in measure_list:
        ms_control_list_dict = dict()
        ms_ipsi_list_dict = dict()
        ms_contra_list_dict = dict()
        for side in side_list:
            # read values for controls, ipsi patients and contra patients
            #-- controls
            ms_control_list_dict[side] = []
            for control in control_list:
                ms_control_list_dict[side].append(
                    measure_dict[measure][control][side])
            #-- ipsi patients
            ms_ipsi_list_dict[side] = []
            for ipsi in patient_list_dict2[side]['ipsi']:
                ms_ipsi_list_dict[side].append(
                    measure_dict[measure][ipsi][side])
            #-- contra patients
            ms_contra_list_dict[side] = []
            for contra in patient_list_dict2[side]['contra']:
                ms_contra_list_dict[side].append(
                    measure_dict[measure][contra][side])
        # generate graph
        fig, ax = plt.subplots(1)
        #-- right hemisphere
        #---- controls
        ax.scatter(
            1*np.ones_like(ms_control_list_dict['right']),
            ms_control_list_dict['right'],
            marker='o',
            c=colour_dict2['right']['control'],
            edgecolors='black',
            s=30,
            label='controls')
        #---- contra-lateral patients
        ax.scatter(
            2*np.ones_like(ms_contra_list_dict['right']),
            ms_contra_list_dict['right'],
            marker='^',
            c=colour_dict2['right']['contra'],
            edgecolors='black',
            s=30,
            label='contralateral')
        #---- ipsi-lateral patients
        ax.scatter(
            3*np.ones_like(ms_ipsi_list_dict['right']),
            ms_ipsi_list_dict['right'],
            marker='d',
            c=colour_dict2['right']['ipsi'],
            edgecolors='black',
            s=30,
            label='ipsilateral')
        #-- left hemisphere
        #---- controls
        ax.scatter(
            4*np.ones_like(ms_control_list_dict['left']),
            ms_control_list_dict['left'],
            marker='o',
            c=colour_dict2['left']['control'],
            edgecolors='black',
            s=30,
            label='controls')
        #---- contra-lateral patients
        ax.scatter(
            5*np.ones_like(ms_contra_list_dict['left']),
            ms_contra_list_dict['left'],
            marker='^',
            c=colour_dict2['left']['contra'],
            edgecolors='black',
            s=30,
            label='contralateral')
        #---- ipsi-lateral patients
        ax.scatter(
            6*np.ones_like(ms_ipsi_list_dict['left']),
            ms_ipsi_list_dict['left'],
            marker='d',
            c=colour_dict2['left']['ipsi'],
            edgecolors='black',
            s=30,
            label='ipsilateral')
        #-- common
        #---- plot limits
        xmin = 0
        xmax = 7
        ymin = 0
        ymax = 0
        if measure == 'volume':
            ymax = 700
        if measure == 'avgthickness':
            if not laplacian:
                ymax = 2.0
            else:
                ymax = 1.2
        ax.axis([xmin, xmax, ymin, ymax])
        #---- labels
        plt_xlabel = ''
        if measure == 'volume':
            plt_xlabel = volume_label
        if measure == 'avgthickness':
            plt_xlabel = avgthick_label
        plt.xlabel(plt_xlabel)
        #---- move volume Y axis to the right
        if measure == 'volume':
            ax.yaxis.set_label_position('right')
            ax.yaxis.tick_right()
        #---- remove tick labels on the X axis
        ax.set_xticklabels([])
        #---- vertical line
        plt.vlines(3.5, ymin, ymax, colors='black', linewidth=1)
        #---- text: left/right hemisphere
        font = dict()
        font['family'] = 'sans-serif'
        font['weight'] = 'bold'
        #------ right
        plt.text(1, ymax/15, 'right hemisphere', fontdict=font)
        #------ left
        plt.text(1+3.5, ymax/15, 'left hemisphere', fontdict=font)
        #---- legend
        control_element = matplotlib.lines.Line2D(
            [0],
            [0],
            marker='o',
            color='white',
            label='controls',
            markerfacecolor='white',
            markeredgecolor='black')
        contra_element = matplotlib.lines.Line2D(
            [0],
            [0],
            marker='^',
            color='white',
            label='contralateral',
            markerfacecolor='white',
            markeredgecolor='black')
        ipsi_element = matplotlib.lines.Line2D(
            [0],
            [0],
            marker='d',
            color='white',
            label='ipsilateral',
            markerfacecolor='white',
            markeredgecolor='black')
        legend_elements = [control_element, contra_element, ipsi_element]
        volume_condition = (measure == 'volume')
        avgthick_condition = (measure == 'avgthickness') and (avgthick_legend)
        if volume_condition or avgthick_condition:
            ax.legend(handles=legend_elements, loc=1)

        # save figure
        fig_path = './scatterplot_{0}.pdf'.format(measure)
        fig_path = os.path.join(
            output_data_path,
            'scatterplot_{0}.pdf'.format(measure))
        plt.savefig(fig_path, bbox_inches='tight')
