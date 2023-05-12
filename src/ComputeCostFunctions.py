import numpy as np
import pandas as pd
import configparser, itertools
from costFunctions import kinematic_cost_function, static_cost_function, weighted_cost_function
from triangulate import triangular_projection
from vtk_utils import read_pvtu


def unpack_str_list(list, dtype=float):
    return [dtype(k) for k in list.split(',')]
    

def compute_kine_cost(result_folder, config):
    Expe_data = config['Experimental Data']
    cost_options = config['Cost Function']
    DIC_data = Expe_data['DIC data']
    
    if config.has_option('Experimental Data', 'DIC time steps'):
        dic_time_steps = unpack_str_list(Expe_data['DIC time steps'], dtype=int)
    else:
        # Read time steps from prm file template
        prm_file = result_folder + ".prm"
        config_template = configparser.ConfigParser()
        with open(prm_file) as fp:
            config_template.read_file(itertools.chain(['[global]'], fp), source=prm_file)
        dt = float(config_template['global']['set Time increments'])
        output_table = unpack_str_list(config_template['global']['set Tabular Time Output Table'])
        dic_time_steps = np.array(output_table)/dt - 1
        dic_time_steps = dic_time_steps.astype('int')

    chi_u = 0
    n_steps = len(dic_time_steps)
    for step, step_simu in enumerate(dic_time_steps):
        pvtu_fname = "{}/solution-{:04}.pvtu".format(result_folder, step_simu)
        nodes, u_SIM = read_pvtu(pvtu_fname)
        if nodes is None:
            # It seems that the simulation has failed, raise penalty value
            chi_u = float(cost_options['penalty'])
            n_steps = 1
            break
        else:
            # Read DIC data
            DIC_step = '{}{}.csv'.format(DIC_data, str(step+1))
            DIC_vals = np.loadtxt(DIC_step)
            pts_DIC = DIC_vals[:, :2]
            u_DIC = DIC_vals[:, 2:4]

            # Associate each DIC measurement to a unique node
            u_SIM_tri, inside_mesh = triangular_projection(nodes, u_SIM, pts_DIC)
            
            # Remove DIC locations outside the RoI
            u_SIM_tri = u_SIM_tri[inside_mesh]

            wgt_str = 'weight by correlation coefficients'
            if config.has_option('Cost Function', wgt_str) and config.getboolean('Cost Function', wgt_str):
                C = DIC_vals[inside_mesh, 4]
            else:
                C = None
            chi_u += kinematic_cost_function(u_SIM_tri, u_DIC[inside_mesh], weights=C)

    return chi_u/n_steps


def compute_stat_cost(result_folder, config):
    # Read data from experimental tensile curve
    tensileCurve = config['Experimental Data']['tensile curve']
    Exper_curve = np.loadtxt(tensileCurve)
    elon_expe = Exper_curve[:, 0]
    stress_expe = Exper_curve[:, 1]

    # Fetch simulated tensile curve
    if config.has_option('Experimental Data', 'tensile direction'):
        tensile_dir = config['Experimental Data']['tensile direction'].lower()
    else:
        tensile_dir = 'x'
    col_E = 'E' + 2 * tensile_dir
    col_sigma = 'T' + 2 * tensile_dir
    try:
        stressstrain = pd.read_csv(result_folder + '/stressstrain.txt', sep='\t')
        gl_eps = stressstrain[col_E]
        elon_simu = -1 + np.sqrt(1 + 2 * gl_eps)  # Compute elongation from Green-Lagrangian strain
        stress_simu = stressstrain[col_sigma]
        return static_cost_function(elon_expe, stress_expe, elon_simu, stress_simu)
    except FileNotFoundError:
        return float(config['Cost Function']['penalty'])


def compute_weighted_cost(result_folder, config):
    # Compute cost functions
    chi_u = compute_kine_cost(result_folder, config)
    chi_F = compute_stat_cost(result_folder, config)
   
    # Return both cost functions, plus the weighted mean
    w_sigma = float(config['Cost Function']['weight on tensile curve'])
    chi = weighted_cost_function(chi_F, chi_u, w1=w_sigma)
    return {'chi_u': chi_u, 'chi_f': chi_F, 'chi': chi}
