# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 21:16:28 2021

@author: Dorian
"""
import configparser
import importlib
import json
import os
import shutil

import numpy as np
import pandas as pd

from CfgGenerator import CfgGenerator
from ComputeCostFunctions import compute_weighted_cost, unpack_str_list


def parse_optional_param(config, section):
    if config.has_section(section):
        options = dict(config[section])
        for key, value in options.items():
            if key in ['forward', 'verbose', 'loginfo', 'time']:
                options[key] = config.getboolean(section, key)
            else:
                options[key] = json.loads(value)
        return options
    else:
        return None


def read_bounds(config):
    lb = np.array(unpack_str_list(config.get("Bounds", "lower")))
    ub = np.array(unpack_str_list(config.get("Bounds", "upper")))
    return lb, ub


def remove_data(path, debug=True):
    if not os.path.isfile(path):
        path = os.path.join(path, '')
    if debug:
        print("delete {}".format(path))
    else:
        if os.path.isfile(path):
            os.remove(path)
        else:
            try:
                shutil.rmtree(path)
            except OSError as e:
                print("Error: %s - %s." % (e.filename, e.strerror))


def optimize(config_file='Config.ini'):
    """
    Run optimization loops.
    
    Parameters
    ----------
    config_file : string, optional
        Path to configuration file. The default is 'Config.ini'.

    Returns
    -------
    scipy.optimize._optimize.OptimizeResult
        Results from optimization.
    """

    # Read optimizer options from config file
    config = configparser.ConfigParser()
    config.read(config_file)

    # Infer the number of investigated parameters from the initial guess
    x0 = np.array([float(i) for i in config['Initial Guess'].values()])

    # Load bounds and check consistency with the initial guess
    lb, ub = read_bounds(config)
    error_msg = 'The size of the {} bounds must be the same as the number of values in initial guess.'
    if len(lb) != len(x0):
        raise ValueError(error_msg.format('upper'))
    if len(ub) != len(x0):
        raise ValueError(error_msg.format('lower'))

    # The minimizer uses absolute differences to compute the gradient. Hence
    # it is a good habit to normalize the parameters so that the investigated 
    # space is an hypercube of size 1.
    x0n = (x0 - lb) / (ub - lb)
    boundsn = np.concatenate((np.zeros((len(x0), 1)), np.ones((len(x0), 1))), axis=1).tolist()

    # Read L-BFGS-B options
    minimize_opt = parse_optional_param(config, 'Minimize')
    kwargs = {'bounds': boundsn, 'args': config, 'options': minimize_opt}
    if config.has_section('Minimize parallel') and config.getboolean('Minimize parallel', 'use parallel minimizer'):
        config['Minimize parallel'].pop('use parallel minimizer')
        kwargs['parallel'] = parse_optional_param(config, 'Minimize parallel')
        module = importlib.import_module('optimparallel')
        minimizer = module.minimize_parallel
    else:
        module = importlib.import_module('scipy.optimize')
        minimizer = module.minimize

    res = minimizer(run_and_compare, x0n, **kwargs)

    # Restore the optimized parameters into un-normalized form
    k = ub - lb
    inv = 1 / k
    res.x = lb + res.x * k
    res.jac = res.jac * inv

    # Plus, provide the estimated hessian
    h = np.linalg.inv(res.hess_inv.todense())
    res.hess = h * np.outer(inv, inv)

    return res


def run_and_compare(pn, config):
    """
    From normalized array of parameters, generate input file for 
    PRISMS-Plasticity, run it, compute the cost function and return the
    cost function to be optimized.

    Parameters
    ----------
    pn : numpy.ndarray
        Normalized array of parameters to be used for PRISMS-Plasticity. This
        array is normalized with respect to the bounds.
    config : configparser.ConfigParser
        Configuration data, parsed by configparser

    Returns
    -------
    float
        Scalarized cost function to be minimized

    """
    # Restore the parameters into their initial (un-normalized) form
    lb, ub = read_bounds(config)
    dom_size = ub - lb
    p = lb + pn * dom_size
    n_p = len(p)

    # Find the absolute tolerance of normalized parameters
    if config.has_option('Minimize', 'eps'):
        eps_jac = float(config['Minimize']['eps'])
    else:
        # It seems that the default value for the absolute tolerance is 1e-8
        eps_jac = 1e-8

    # First, look in log file if this simulation has been run before
    logfile = config['Log File']['file path']
    chi_names = ['chi_u', 'chi_f', 'chi']
    col_names = list(config['Initial Guess'].keys()) + chi_names
    try:
        # Try to read log file
        prev = pd.read_csv(logfile)
        if prev.size == 0:
            # Only the header is present
            cost = np.array([])
        else:
            mat = prev.to_numpy()
            matn = (mat[:, :n_p] - lb) / dom_size
            existing = np.all(np.isclose(matn, pn, atol=eps_jac / 10, rtol=0.), axis=1)
            cost = mat[existing, -1]
    except (pd.errors.EmptyDataError, FileNotFoundError):
        cost = np.array([])

    if cost.size != 0:
        # If the simulation was run before, just return the related cost function
        return cost[0]
    else:
        # Otherwise, run this simulation and compute the cost function
        # Generate a dictionary from the parameters
        d = dict(zip(config['Initial Guess'].keys(), p))
        prm_name, lh_name, fname = CfgGenerator(d, config)

        # Run simulation and wait till the end
        if config.has_option('Slurm', 'use Slurm') and config.getboolean('Slurm', 'use Slurm'):
            # Use SLURM workload manager
            batch_file = config['Slurm']['batch file']
            cmd = "sbatch --wait {} {}".format(batch_file, prm_name)
        else:
            # Run PRISMS directly
            bin_path = config['PRISMS']['prisms command line']
            cmd = "{} {}".format(bin_path, prm_name)

        if config.has_option('Debug', 'fake simulations') and config.getboolean('Debug', 'fake simulations'):
            execute = print
        else:
            execute = os.system
        execute(cmd)

        # Compute the cost function
        chis = compute_weighted_cost(fname, config)

        # Remove conf files and results
        debug = config.has_option('Debug', 'fake deletions') and config.getboolean('Debug', 'fake deletions')
        remove_data(lh_name, debug=debug)
        remove_data(prm_name, debug=debug)
        remove_data(fname, debug=debug)

        # Append the cost functions to the log file
        a = np.concatenate((p, np.array([chis['chi_u'], chis['chi_f'], chis['chi']])))
        df = pd.DataFrame(data=[a], columns=col_names)
        try:
            _ = pd.read_csv(logfile)
            df.to_csv(logfile, index=False, header=False, mode='a')
        except (pd.errors.EmptyDataError, FileNotFoundError):
            df.to_csv(logfile, index=False)

        # Return the function to be minimized
        return chis['chi']
