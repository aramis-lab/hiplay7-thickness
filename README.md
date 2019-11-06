# 7T Hippocampus Thickness from Diffeomorphic Vector Field

Code to extract central surface representations along with thickness
maps from MR images of the hippocampus at 7 Tesla.

**Authors**: Émilie Gerardin (Brain and Spine Institute), Ana B.
Graciano Fouquier(Brain and Spine Institute), Alexis Guyot (Brain and
Spine Institute).


## Overview

This folder contains all the code to run the experiments presented in
the following article:
> A Guyot\*, AB Graciano Fouquier\*, É Gerardin\*, M Chupin, JA Glaunès,
> L Marrakchi-Kacem, J Germain, C Boutet, C Cury, L Hertz-Pannier,
> A Vignaud, S Durrleman, TR Henry, PF Van de Moortele, A Trouvé,
> O Colliot.
> (\* denotes equal first authors)
> 'A Diffeomorphic Vector Field Approach to Analyze the Thickness of the
> Hippocampus from 7T MRI'



## Installation

The following software and libraries:
- Matlab
- Python (version 3.2 or newer)
    - numpy
    - scipy
    - matplotlib
    - joblib
    - nibabel
    - nilearn
    - vtk
- Deformetrica (version 4.2.0)

In case you are not sure you already have the relevant Python libraries
(numpy, scipy, matplotlib, joblib, nibabel, nilearn, vtk), we recommend
installing Miniconda, a program that lets you install and run Python
packages and their dependencies into local, user-defined environments.

Miniconda can be obtained at the following website:
https://docs.conda.io/en/latest/miniconda.html. 
Please make sure you choose the version corresponding to your operating
system (Windows, Mac OS X or Linux) and to the architecture of your
computer (32bit or 64bit).

Deformetrica (http://www.deformetrica.org/) can be installed on a Conda
environment.
See https://gitlab.com/icm-institute/aramislab/deformetrica for more
information.


## Usage

Extract the archive 'inputData.tar.gz' present at the root of the git
repository with the following command:
```
tar -xvzf inputData.tar.gz
```
This will create a folder 'inputData' at the root of the current
directory. You can move the 'inputData' folder wherever you wish, as
long as you keep track of the location (which will be referred to as
[input_folder] in the following instructions).

**Optional**: if you have installed dependencies via miniconda, activate
the conda environment that contains the dependencies.

To get the results presented in the article, run the following commands:

- Experiment III. A. _Robustness with respect to kernel size and
anisotropy_: 
```
python ./launch_experiment_A.py [input_folder] [output_folder_A]
```
- Experiment III. B. _Influence of inter-rater variability_:
```
python ./launch_experiment_B.py [input_folder] [output_folder_B] (-n [n_cores])
```
- Experiment III. C. _Application to in vivo 7T MRI group studies_:
```
python ./launch_experiment_C.py [input_folder] [output_folder_C] (-n [n_cores])
```

Where:
- [input_folder] is the path to the extracted archive containing the
input data for all the experiments described in the article.
- [output_folder_A], [output_folder_B] and [output_folder_C] are
separate empty folders where the output data for experiments A, B and C
respectively, will be stored.
- n_cores is an optional argument defining the number of cores to be
used in parallel tasks. Default is 1.

## Output

### Experiment III. A. _Robustness with respect to kernel size and anisotropy_

[output_folder_A] (as described in section **Usage**) will be
populated with the following folders:

- '4-centralSurfaceAndThicknessEstimation': contains files for the
central surface calculated with a kernel size t=10

- '5-influenceOfKernelSize': contains two subfolders
    - '1-computeMaps': contains files for the central surfaces computed
    for kernel sizes 3, 5, 10 and 15
    - '2-compareMaps': contains comparisons from the surfaces at kernel
    size 3, 5, 15 to the surface at kernel size 10 using 3 metrics
    (thickness correlation, mean abs. thickness difference and mean
    inter-surface distance).

- '6-influenceOfAnisotropy': contains three subfolders
    - '1-createAnisotropicSegmentations': populated with subfolders
    'subample_1', ..., 'subsample_6', each containing subsampled
    segmentations for subsampling factors in 1, 2, ..., 6
    - '2-computeMaps': contains files for the central surfaces computed
    for each subsampled segmentation
    - '3-compareMaps': contains comparisons from the central surfaces
    obtained with subsampling factors 2, 3, ..., 6 to the central
    surface obtained with subsampling factor 1

### Experiment III. B. _Influence of inter-rater variability_

[output_folder_B] (as described in section **Usage**) will be
populated with the following folders:

- '4-centralSurfaceAndThicknessEstimation': contains files for the
central surface calculated with a kernel size t=10 for rater1/rater2,
subject1/subject2/subject3/subject4 and left/right hippocampus

- '5-interRaterAnalysis': contains two subfolders
    - '1-compareIndivididual': contains the comparison (using the 3
    metrics described above) between rater 1 and rater 2 for the
    left/right hippocampus of subject1/subject2/subject3/subject4
    - '2-compareAverage': contains the three metrics averaged across all
    four subjects for both left and right hippocampus

### Experiment III. C. _Application to in vivo 7T MRI group studies_

[output_folder_C] (as described in section **Usage**) will be
populated with the following folders:

- '4-groupStudy': single folder containing five subfolders
    - '1-centralSurfaceAndThicknessEstimation': contains files for the
    central surfaces calculated with a kernel size of all subjects
    (controls: control1, control2, control3, control4, control6,
    control8, control9, control10, control11 ; left contralateral
    patients: patient1, patient2, patient6 ; left ipsilateral patients:
    patient3, patient5, patient8, patient9, patient12), for both left
    and right hippocampus
    - '2-centralSurfaceTemplates': contain the template hippocampus
    surface computed by Deformetrica for both {controls + contralateral
    patients} and {controls + ipsilateral patients} groups for both left
    and right hippocampus
    - '3-templateProjections': for the two groups described above and
    for both hippocampus, contains the average of all projected
    thickness maps from controls, the average of all projected thickness
    maps from ipsi- or contra- lateral patients and the difference map
    between these two averages
    - '4-avgThicknessVolumeComputation': contains files reporting the
    volume and average thickness of all hippocampi described in this
    experiment
    - '5-spearmanCorrelationsComputation': contains file
    'spearman_correlations.json' which reports Spearman's rank
    correlation coefficents between volume and average thickness for the
    left and right hippocampi (computed across all subjects)

## License

With the exceptions of files taken from external libraries, this code is
released under the terms of the MIT License.
See file LICENSE.txt for further precisions.
