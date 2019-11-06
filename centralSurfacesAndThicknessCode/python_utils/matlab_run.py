#! /usr/bin/python

"""Matlab utils functions

Wraps subprocess.run('matlab ...') into a higher level collection of
functions to allow Matlab code to run within a Python script.
"""

import subprocess


def create_matlab_script(matlab_func, matlab_arg_list, matlab_path_list):
    """Create the Matlab script to be run in interpreter

    Build the matlab script to be run in the Matlab interpreter from
    function name, arguments passed to functions and required paths.

    Args:
        matlab_func (str): name of the function to be run. File
            [matlab_func].m has to belong to the matlabpath (i.e.,
            system function or a file within matlab_path_list)
        matlab_arg_list (list): List of strings. Arguments given to
            matlab_func. The order has to match that of the function
            arguments.
        matlab_path_list (list): List of the paths that must be provided
            to the Matlab interpreter prior to running the function

    Returns
        matlab_script (str): matlab script to be run
            (e.g., 'my_function arg1 arg2')
    """
    # initialise matlab command
    matlab_script = ''

    # Add all paths to the command
    for matlab_path in matlab_path_list:
        matlab_script = '{0}addpath(genpath(\'{1}\')) ; '.format(
            matlab_script, matlab_path)

    # Add Matlab function to the command
    matlab_script = '{0}{1} '.format(matlab_script, matlab_func)

    # Add all the function arguments to the command
    for matlab_arg in matlab_arg_list:
        matlab_script = '{0} {1}'.format(matlab_script, matlab_arg)

    return matlab_script


def matlabscript_to_shellcmd(matlab_script):
    """Create the shell command to run a Matlab script

    Generate a shell command that will launch a Matlab interpreter
    tor run a matlab script.

    Args:
        matlab_script (str): matlab script to be run
            (e.g., 'my_function arg1 arg2') in a Matlab interpreter

    Returns:
        shell_cmd (str): shell command to run matlab script in a
            Matlab interpreter
    """
    # generate shell command to run the script within Matlab
    #-- call to Matlab executable
    matlab_executable = None
    try:
        check_wich_matlab = subprocess.check_output(['which', 'matlab'])
        matlab_executable = check_wich_matlab.decode().rstrip()
    except subprocess.CalledProcessError:
        print('Matlab not found')
    shell_cmd = matlab_executable
    #-- flags
    matlab_cli_flags = ''
    shell_cmd = '{0} {1}'.format(shell_cmd, matlab_cli_flags)
    #---- flags to run script
    run_script_flag = '-r'
    shell_cmd = '{0} {1}'.format(shell_cmd, run_script_flag)
    #-- script called by matlab + exit Matlab interpreter
    shell_cmd = '{0} "{1} ; exit;"'.format(shell_cmd, matlab_script)

    return shell_cmd


def run_matlab_func(matlab_func, matlab_arg_list, matlab_path_list):
    """Run matlab function on arguments

    Build the matlab script to be run in the Matlab interpreter from
    function name, arguments passed to functions and required paths,
    then build shell command to run the script in Matlab interpreter,
    then run the shell command.

    Args:
        matlab_func (str): name of the function to be run. File
            [matlab_func].m has to belong to the matlabpath (i.e.,
            system function or a file within matlab_path_list)
        matlab_arg_list (list): List of strings. Arguments given to
            matlab_func. The order has to match that of the function
            arguments.
        matlab_path_list (list): List of the paths that must be provided
            to the Matlab interpreter prior to running the function

    Returns
        N/A
    """
    # build matlab script
    matlab_script = create_matlab_script(
        matlab_func, matlab_arg_list, matlab_path_list)

    # build shell command to run matlab script in Matlab interpreter
    shell_cmd = matlabscript_to_shellcmd(matlab_script)

    # run the shell command
    print('run')
    print(shell_cmd)
    subprocess.run(shell_cmd, shell=True)
