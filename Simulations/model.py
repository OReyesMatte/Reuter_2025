import numpy as np
import multiprocessing as mp
from ddeint import ddeint

# default model parameters
model_parameters = {
    'phi_WT': 4*10**-9,      # adsorption rate ancestor (min^-1)
    'phi_LP': 2*10**-9,      # adsorption rate large-plaque mutant (min^-1)
    'phi_SP': 6*10**-9,      # adsorption rate small-plaque mutant (min^-1)
    'beta_WT': 200,          # burst size ancestor (PFU)
    'beta_LP': 200,          # burst size large-plaque mutant (PFU)
    'beta_SP': 200,          # burst size small-plaque mutant (PFU)
    'tau': 17,               # lysis time (min)
    'd_WT': 1.7*10**-2,      # decay rate ancestor (min^-1)
    'd_LP': 1.4*10**-2,      # decay rate large-plaque mutant (min^-1)
    'd_SP': 6*10**-3,        # decay rate small-plaque mutant (min^-1)
    'm': 0,                  # mutation rate, if using mutations (min^-1)
    'initial_host': 10**8,   # initial host density
    'initial_MOI': 0.1,      # initial MOI (initial phage density = initial_MOI * initial_host)
    'transfer_volume': 0.02, # fraction of the volume that is transferred
    'transfer_times': 5*[30] # sequence of transfer times
}

# parameters for the numerical solution
numerical_parameters = {
    'dt': 10**-2
}

# the model of bacteria-phage interaction
def bacteria_phage_model(Y, t, p):  
    
    # densities at time t
    B, Bi, P_WT, P_LP, P_SP = Y(t)
    
    # densities at time t - tau in the past (delay)
    B_t, Bi_t, P_WT_t, P_LP_t, P_SP_t = Y(t - p['tau'])
    
    # bacteria
    dBdt = - (p['phi_WT']*P_WT + p['phi_LP']*P_LP + p['phi_SP']*P_SP)*B
    
    # lysis is determined by adsorptions tau minutes ago (delayed state)
    lysis_WT = B_t*P_WT_t
    lysis_LP = B_t*P_LP_t                                                                                                     
    lysis_SP = B_t*P_SP_t
    
    # if within the first tau minutes of transfer, take dilution into account for delayed state
    if t < p['tau']:
        lysis_WT = p['transfer_volume']*lysis_WT
        lysis_LP = p['transfer_volume']*lysis_LP
        lysis_SP = p['transfer_volume']*lysis_SP
        
    # infected bacteria
    dBidt = p['phi_WT']*(B*P_WT - lysis_WT) + p['phi_LP']*(B*P_LP - lysis_LP) + p['phi_SP']*(B*P_SP - lysis_SP)
    
    # three phage types (wild type/ancestor WT, large-plaque mutant LP, small-plaque mutant SP)
    dP_WTdt = p['phi_WT']*((1-p['m'])*p['beta_WT']*lysis_WT - (B+Bi)*P_WT) + 0.5*p['m']*p['beta_LP']*p['phi_LP']*lysis_LP + 0.5*p['m']*p['beta_SP']*p['phi_SP']*lysis_SP - p['d_WT']*P_WT
    dP_LPdt = p['phi_LP']*((1-p['m'])*p['beta_LP']*lysis_LP - (B+Bi)*P_LP) + 0.5*p['m']*p['beta_WT']*p['phi_WT']*lysis_WT + 0.5*p['m']*p['beta_SP']*p['phi_SP']*lysis_SP - p['d_LP']*P_LP
    dP_SPdt = p['phi_SP']*((1-p['m'])*p['beta_SP']*lysis_SP - (B+Bi)*P_SP) + 0.5*p['m']*p['beta_WT']*p['phi_WT']*lysis_WT + 0.5*p['m']*p['beta_LP']*p['phi_LP']*lysis_LP - p['d_SP']*P_SP
    
    return np.array([dBdt, dBidt, dP_WTdt, dP_LPdt, dP_SPdt])

# solves a single transfer experiment with a given transfer protocol and set of parameters
def run_experiment(transfer_experiment):
    print(f"process {(mp.current_process().pid)} running experiment " + transfer_experiment['experiment_id'] + '\n')

    # initial condition and history function for the first transfer period
    transfer_initial_condition = lambda t: [0]*len(transfer_experiment['initial_condition']) if t<0 else transfer_experiment['initial_condition']

    solution = np.empty((0, 1+len(transfer_experiment['initial_condition'])))
    for transfer_i, t_end in enumerate(transfer_experiment['parameters']['transfer_times']):
        
        # the time interval for which to solve the dde [0, transfer_time)
        time = np.linspace(0, t_end, int(t_end/numerical_parameters['dt']))
        
        # call ddeint to numerically solve the dde governing the continuous within-transfer dynamics
        transfer_solution = ddeint(bacteria_phage_model, transfer_initial_condition, time, fargs=(transfer_experiment['parameters'],))
        
        # define the initial condition for the next transfer (fresh bacterial host + the last values from the previous transfer)
        transfer_initial_condition = lambda t: transfer_solution[int(t/numerical_parameters['dt'])] if t<0 else np.hstack([[transfer_experiment['parameters']['transfer_volume']*transfer_solution[:, 0][-1] + transfer_experiment['initial_condition'][0]], transfer_experiment['parameters']['transfer_volume']*np.array([transfer_solution[:, i][-1] for i in range(1, len(transfer_experiment['initial_condition']))])])
        
        # after every 180 min transfer inoculate next transfer with MOI 0.1
        if transfer_experiment['parameters']['transfer_times'][transfer_i] == 180:

            final_WT = transfer_solution[:, 2][-1]
            final_LP = transfer_solution[:, 3][-1]
            final_SP = transfer_solution[:, 4][-1]
            
            final_P = final_WT + final_LP + final_SP
            
            init_WT = (final_WT/final_P)*transfer_experiment['parameters']['initial_MOI']*transfer_experiment['initial_condition'][0]
            init_LP = (final_LP/final_P)*transfer_experiment['parameters']['initial_MOI']*transfer_experiment['initial_condition'][0]
            init_SP = (final_SP/final_P)*transfer_experiment['parameters']['initial_MOI']*transfer_experiment['initial_condition'][0]
            
            transfer_initial_condition = lambda t: [0]*len(transfer_experiment['initial_condition']) if t<0 else [transfer_experiment['initial_condition'][0], 0, init_WT, init_LP, init_SP]
        
        # add the current transfer period
        solution = np.vstack([solution, np.c_[time, transfer_solution]])

    return solution

# start multiple parallel processes, one for each defined transfer experiment
def run_parallel_experiments(transfer_experiments):
    solutions = []
    p = mp.Pool()
    for result in p.map(run_experiment, transfer_experiments):
        solutions.append(result)
    
    return solutions
