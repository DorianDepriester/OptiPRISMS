import numpy as np
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
    
    if config.has_option('Experimental Data', 'time steps'):
        time_steps = unpack_str_list(Expe_data['time steps'], dtype=int)
    else:
        # Read time steps from prm file template
        prm_file = result_folder + ".prm"
        config_template = configparser.ConfigParser()
        with open(prm_file) as fp:
            config_template.read_file(itertools.chain(['[global]'], fp), source=prm_file)
        dt = float(config_template['global']['set Time increments'])
        output_table = unpack_str_list(config_template['global']['set Tabular Time Output Table'])
        time_steps = np.array(output_table)/dt - 1
        time_steps = time_steps.astype('int')

    # Read DIC data and check consistency with simulation parameters
    dataDIC = np.loadtxt(DIC_data)
    pts_DIC = dataDIC[:,:2]
    good_nbr = (2 + 3*len(time_steps))
    if dataDIC.shape[1] != good_nbr:
        err_msg = 'The number of time steps ({}) is inconsistent with the number of columns in {} ({}, instead of {})'\
            .format(len(time_steps), DIC_data, dataDIC.shape[1], good_nbr)
        raise ValueError(err_msg)

    chi_u = 0
    n_steps = len(time_steps)
    for step, step_simu in enumerate(time_steps):
        pvtu_fname = "{}/solution-{:04}.pvtu".format(result_folder, step_simu)
        nodes, u_SIM = read_pvtu(pvtu_fname)
    
        if nodes is None:
            # It seems that the simulation has failed, raise penalty value
            chi_u = float(cost_options['penalty'])
            n_steps = 1
            break
            
        else:
            # Associate each DIC measurement to a unique node
            u_SIM_tri, inside_mesh = triangular_projection(nodes, u_SIM, pts_DIC)
            
            # Remove DIC locations which are too far from mesh nodes
            u_SIM_tri = u_SIM_tri[inside_mesh]
            u_DIC = dataDIC[inside_mesh, (3*step+2):(3*step+4)]
            C = dataDIC[inside_mesh, 3*step+4]
            chi_u += kinematic_cost_function(u_SIM_tri[:, :2], u_DIC, weights=C)
  
    return chi_u/n_steps


def compute_stat_cost(result_folder, config):
    tensileCurve = config['Experimental Data']['tensile curve']
    try:
        stressstrain = np.loadtxt(result_folder + '/stressstrain.txt', skiprows=1)
        
        # PRISMS-Plasticity returns the Green-Lagrangian strain
        Exx = stressstrain[:, 0]
        eps_xx_simu = -1+np.sqrt(1+2*Exx)   # Compute elongation from Exx
        Sxx_simu = stressstrain[:, 6]
    
        # Read data from experimental tensile curve
        Exper_curve = np.loadtxt(tensileCurve)
        eps_xx_expe = Exper_curve[:, 0]
        Sxx_expe = Exper_curve[:, 1]
        return static_cost_function(eps_xx_expe, Sxx_expe, eps_xx_simu, Sxx_simu)
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
