import numpy as np
from ddeint import ddeint
import multiprocessing as mp
from matplotlib import pyplot as pp

# main model parameters
main_parameters = {
    'phi_WT': 4.15*10**-9,  # adsorption rate ancestor (min^-1)
    'phi_LP': 9.87*10**-10, # adsorption rate large-plaque mutant (min^-1)
    'phi_SP': 6*10**-9,     # adsorption rate small-plaque mutant (min^-1)
    'beta_WT': 200,         # burst size ancestor (PFU)
    'beta_LP': 200,         # burst size large-plaque mutant (PFU)
    'beta_SP': 200,         # burst size small-plaque mutant (PFU)
    'tau': 17,              # lysis time (min)
    'd_WT': 1.7*10**-2,     # decay rate ancestor (min^-1)
    'd_LP': 1.4*10**-2,     # decay rate large-plaque mutant (min^-1)
    'd_SP': 6*10**-3,       # decay rate small-plaque mutant (min^-1)
    'm': 10**-5             # mutation rate (min^-1)
}

# the model of bacteria-phage interaction
def bacteria_phage_model(Y, t, p):  
    # densities at current time t
    B, Bi, P_WT, P_LP, P_SP = Y(t)
    
    # densities at time t - tau in the past (delay)
    B_t, Bi_t, P_WT_t, P_LP_t, P_SP_t = Y(t - p['tau'])
    
    # dynamics of un-infected bacteria
    dBdt = - (p['phi_WT']*P_WT + p['phi_LP']*P_LP + p['phi_SP']*P_SP)*B
    
    # dynamics of infected bacteria
    dBidt = (p['phi_WT']*(B*P_WT - B_t*P_WT_t) + p['phi_LP']*(B*P_LP - B_t*P_LP_t) + p['phi_SP']*(B*P_SP - B_t*P_SP_t))
    
    # dynamics of the three phage types (wild type/ancestor WT, large-plaque mutant LP, small-plaque mutant SP)
    dP_WTdt = p['phi_WT']*((1-p['m'])*p['beta_WT']*B_t*P_WT_t - (B+Bi)*P_WT) + 0.5*p['m']*p['beta_LP']*p['phi_LP']*B_t*P_LP_t + 0.5*p['m']*p['beta_SP']*p['phi_SP']*B_t*P_SP_t - p['d_WT']*P_WT
    dP_LPdt = p['phi_LP']*((1-p['m'])*p['beta_LP']*B_t*P_LP_t - (B+Bi)*P_LP) + 0.5*p['m']*p['beta_WT']*p['phi_WT']*B_t*P_WT_t + 0.5*p['m']*p['beta_SP']*p['phi_SP']*B_t*P_SP_t - p['d_LP']*P_LP
    dP_SPdt = p['phi_SP']*((1-p['m'])*p['beta_SP']*B_t*P_SP_t - (B+Bi)*P_SP) + 0.5*p['m']*p['beta_WT']*p['phi_WT']*B_t*P_WT_t + 0.5*p['m']*p['beta_LP']*p['phi_LP']*B_t*P_LP_t - p['d_SP']*P_SP
    
    return np.array([dBdt, dBidt, dP_WTdt, dP_LPdt, dP_SPdt])


# solves a single transfer experiment with a given transfer protocol and set of parameters
def run_experiment(transfer_experiment):
    print(f"process {(mp.current_process().pid)} running transfer experiment " + transfer_experiment['experiment_id'] + '\n')

    # initial condition and history function for the first transfer period
    init_transfer = lambda t: [0]*len(transfer_experiment['initial_condition']) if t<0 else transfer_experiment['initial_condition']
    
    dt = 10**-2
    all_transfers = np.empty((0,6))
    for t_end in transfer_experiment['transfer_times']:
        
        # the time interval for which to solve the dde [0, transfer_time)
        time = np.linspace(0, t_end, int(t_end/dt))
        
        # call ddeint to numerically solve the dde governing the continuous within-transfer dynamics
        transfer = ddeint(bacteria_phage_model, init_transfer, time, fargs=(transfer_experiment['parameters'],))
        
        # define the initial condition (fresh bacterial host + the last values from the previous transfer)
        # and history function for the next transfer period (discrete step/map)
        init_transfer = lambda t: transfer_experiment['transfer_vol']*transfer[int(t/dt)] if t<0 else np.hstack([[transfer_experiment['transfer_vol']*transfer[:, 0][-1] + transfer_experiment['initial_condition'][0]], transfer_experiment['transfer_vol']*np.array([transfer[:, 1][-1], transfer[:, 2][-1], transfer[:, 3][-1], transfer[:, 4][-1]])])
        
        # save the current transfer period
        all_transfers = np.vstack([all_transfers, np.c_[time, transfer]])

    return all_transfers

# start multiple parallel processes, one for each defined transfer experiment
def run_parallel_experiments(transfer_experiments):
    solutions = []
    p = mp.Pool()
    for result in p.map(run_experiment, transfer_experiments):
        solutions.append(result)
    
    return solutions


def plot_transfer_experiments(transfer_experiments, solutions):
    for i, solution in enumerate(solutions):
        sol_t = solution[:, 0]        
        sol_WT = solution[:, 3]
        sol_LP = solution[:, 4]
        sol_SP = solution[:, 5]
        
        transfer_i = np.where(sol_t == 0)[0]
        
        fig, ax = pp.subplots(1, 1, sharex=False, figsize=(5, 3))
        
        ax.set_xlabel('Transfer')
        ax.set_ylabel('Genotype frequency')
        
        ticks_label = [1] + [ti+1 if (ti+1)%5==0 else '' for ti in range(1, len(transfer_i))]    
            
        ax.set_xticks(transfer_i)
        ax.set_xticklabels(ticks_label)
        
        ax.plot(transfer_i, sol_WT[transfer_i]/(sol_WT[transfer_i] + sol_LP[transfer_i] + sol_SP[transfer_i]), 'o--', lw=2, color='#29485d', label='ancestor')
        ax.plot(transfer_i, sol_LP[transfer_i]/(sol_WT[transfer_i] + sol_LP[transfer_i] + sol_SP[transfer_i]), 'o--', lw=2, color='#ae017e', label='large-plaque mutant')
        ax.plot(transfer_i, sol_SP[transfer_i]/(sol_WT[transfer_i] + sol_LP[transfer_i] + sol_SP[transfer_i]), 'o--', lw=2, color='#c6d325', label='small-plaque mutant')
    
        pp.tight_layout()
        pp.savefig(transfer_experiments[i]['experiment_id'] + ".pdf", bbox_inches='tight')

def plot_optimality_figure(transfer_experiments, solutions):
    
    results = []
    for i, result in enumerate(solutions):
        results.append([transfer_experiments[i]['parameters']['phi_LP'], solutions[i]])
    
    pp.figure()
    colors = ['#29485d', '#ae017e', '#c6d325']
    pp.plot([vv[0] for vv in results], [vv[1][-1][4]/vv[1][-1][3] for vv in results], lw=3, color='black')
    pp.axhline(1, color='grey', lw=2, alpha=0.5)
    pp.vlines([main_parameters['phi_WT'], main_parameters['phi_LP']], -1, [1, 13.2], color=[colors[0], colors[1]], zorder=1)
    pp.yticks(list(pp.yticks()[0]) + [1])
    pp.xlabel(r'Mutant adsorption rate (min$^{-1}$)')
    pp.ylabel(r'Relative free phage particles after 30 mins')
    pp.text(1.1*main_parameters['phi_WT'], 1.5, 'Ancestor', color=colors[0], ha='left')
    pp.text(1.2*main_parameters['phi_LP'], 13.5, 'Large plaque mutant', color=colors[1], ha='left')
    pp.xscale('log')
    pp.ylim(-1, 15)
    
    pp.tight_layout()
    pp.savefig("Fig_4B.pdf", bbox_inches='tight')


# creates Figs. 4A, 5C and Suppl. Fig. S5
def make_transfer_figures():
    
    # parameter set without decay (Figs. 4A and S5A)
    parameters_no_decay = main_parameters.copy()
    parameters_no_decay['d_WT'] = 0
    parameters_no_decay['d_LP'] = 0
    parameters_no_decay['d_SP'] = 0
    
    # parameter set with equal adsorption rates (Fig. S5B)
    parameters_same_phi = main_parameters.copy()
    parameters_same_phi['d_LP'] = parameters_same_phi['phi_WT']
    parameters_same_phi['d_SP'] = parameters_same_phi['phi_WT']

    # set initial condition for first transfer
    MOI = 0.1
    initial_host = 10**8
    initial_phage = MOI*initial_host
    initial_conditions = [initial_host, 0, initial_phage, 0, 0]
    
    # transfer volume (relative to total volume)
    transfer_vol = 0.02
    
    # number of transfers
    n_transfers = 50
    
    # the durations of each transfer for the different transfer experiments (in minutes)
    # 30mins: [30, 30, 30, 30, 180, ...]
    transfer_times_30mins = [180 if (i+1)%5==0 else 30 for i in range(n_transfers)]
    # 3hrs: [180, 180, 180, ...]
    transfer_times_3hrs = n_transfers*[180]
    
    # this defines the specific transfer experiments we want to solve
    # each experiment consists of a set of biological parameters (adsorption rate etc.),
    # the transfer protocol (times & transferred fraction) and the initial condition for the very first transfer period
    transfer_experiments = list([
                                {
                                    'experiment_id': 'Fig_5C_1',
                                    'parameters': main_parameters,
                                    'transfer_times': transfer_times_30mins,
                                    'transfer_vol': transfer_vol,
                                    'initial_condition': initial_conditions
                                },
                                {
                                    'experiment_id': 'Fig_5C_2',
                                    'parameters': main_parameters,
                                    'transfer_times': transfer_times_3hrs,
                                    'transfer_vol': transfer_vol,
                                    'initial_condition': initial_conditions
                                },
                                {
                                    'experiment_id': 'Fig_4A',
                                    'parameters': parameters_no_decay,
                                    'transfer_times': transfer_times_30mins,
                                    'transfer_vol': transfer_vol,
                                    'initial_condition': initial_conditions
                                },
                                {
                                    'experiment_id': 'Fig_S5A',
                                    'parameters': parameters_no_decay,
                                    'transfer_times': transfer_times_3hrs,
                                    'transfer_vol': transfer_vol,
                                    'initial_condition': initial_conditions
                                },
                                {
                                    'experiment_id': 'Fig_S5B',
                                    'parameters': parameters_same_phi,
                                    'transfer_times': transfer_times_3hrs,
                                    'transfer_vol': transfer_vol,
                                    'initial_condition': initial_conditions
                                }
                            ])
    
    # run all defined transfer experiments in parallel
    solutions = run_parallel_experiments(transfer_experiments)
    
    # plot the results from each transfer experiment
    plot_transfer_experiments(transfer_experiments, solutions)


# creates Fig. 4B
def make_figure_4B():

    # the range of mutant adsorption rates competing against the ancestor
    phi_range = np.logspace(-11.3, -8, base=10, num=100)
    phi_range = np.sort(np.append(phi_range, [main_parameters['phi_WT'], main_parameters['phi_LP']]))

    # initial condition
    MOI = 0.1
    initial_host = 10**8
    initial_phage = MOI*initial_host
    initial_conditions = [initial_host, 0, 0.5*initial_phage, 0.5*initial_phage, 0]
  
    # define a transfer experiment for each mutant adsorption rate
    transfer_experiments = []
    for i, phi_LP in enumerate(phi_range):
        parameters_var_phi = main_parameters.copy()
        parameters_var_phi['phi_LP'] = phi_LP
        transfer_experiments.append(
                                    {
                                        'experiment_id': 'Fig_4B_' + str(i+1),
                                        'parameters': parameters_var_phi,
                                        'transfer_times': [30],
                                        'transfer_vol': 0.02,
                                        'initial_condition': initial_conditions
                                    }
                                )
    # run all defined transfer experiments in parallel
    solutions = run_parallel_experiments(transfer_experiments)
    
    # plot the results (Fig. 4B)
    plot_optimality_figure(transfer_experiments, solutions)
    

make_transfer_figures()
make_figure_4B()
