import numpy as np
from matplotlib import pyplot as pp
from model import model_parameters, run_parallel_experiments


def make_figure_4C():
    
    n_transfers = 20
    transfer_times = 180
    
    # parameters for figure 4C_1
    parameters_adsorption_only = model_parameters.copy()
    parameters_adsorption_only['transfer_times'] = n_transfers*[transfer_times]
    parameters_adsorption_only['d_WT'] = 0
    parameters_adsorption_only['d_LP'] = 0
    parameters_adsorption_only['d_SP'] = 0
    
    # parameters for figure 4C_2
    parameters_decay_only = model_parameters.copy()
    parameters_decay_only['transfer_times'] = n_transfers*[transfer_times]
    parameters_decay_only['phi_LP'] = parameters_decay_only['phi_WT']
    parameters_decay_only['phi_SP'] = parameters_decay_only['phi_WT']
    
    # initial condition for first transfer
    initial_host = model_parameters['initial_host']
    initial_phage = model_parameters['initial_MOI']*initial_host
    initial_condition = [initial_host, 0, 0.8*initial_phage, 0.1*initial_phage, 0.1*initial_phage]
    
    # this defines the specific transfer experiments we want to solve
    # each experiment consists of a set of biological/technical parameters (adsorption rate, burst size, transfer protocol, etc.),
    # and the initial condition for the very first transfer period
    transfer_experiments = list([
                                {
                                    'experiment_id': 'Fig_4C_1',
                                    'parameters': parameters_adsorption_only,
                                    'initial_condition': initial_condition
                                },
                                {
                                    'experiment_id': 'Fig_4C_2',
                                    'parameters': parameters_decay_only,
                                    'initial_condition': initial_condition
                                }
                            ])
    
    # run all defined transfer experiments in parallel
    solutions = run_parallel_experiments(transfer_experiments)
    
    # plot the results for each experiment
    for i, solution in enumerate(solutions):
        sol_t = solution[:, 0]        
        sol_WT = solution[:, 3]
        sol_LP = solution[:, 4]
        sol_SP = solution[:, 5]
        
        transfer_i = np.where(sol_t == 0)[0]
        transfer_i = np.append((transfer_i-1)[1:], len(sol_t)-1)
        
        fig, ax = pp.subplots(1, 1, sharex=False, figsize=(5, 3))
        
        ax.set_xlabel('Transfer')
        ax.set_ylabel('Genotype frequency')
        
        x_tick_labels = [ti+1 if (ti+1)%5==0 else '' for ti in range(0, len(transfer_i))]
        ax.set_xticks(transfer_i)
        ax.set_xticklabels(x_tick_labels)
        ax.set_ylim(-0.05, 1.05)

        ax.plot(transfer_i, sol_LP[transfer_i]/(sol_WT[transfer_i] + sol_LP[transfer_i] + sol_SP[transfer_i]), 'o--', lw=2, color='#ae017e', label=r'slow adsorption ($\phi = 2 \cdot 10^{-9}$)')
        ax.plot(transfer_i, sol_WT[transfer_i]/(sol_WT[transfer_i] + sol_LP[transfer_i] + sol_SP[transfer_i]), 'o--', lw=2, color='#29485d', label=r'intermediate adsorption ($\phi = 4 \cdot 10^{-9}$)')
        ax.plot(transfer_i, sol_SP[transfer_i]/(sol_WT[transfer_i] + sol_LP[transfer_i] + sol_SP[transfer_i]), 'o--', lw=2, color='#c6d325', label=r'fast adsorption ($\phi = 6 \cdot 10^{-9}$)')
        
        pp.tight_layout()
        pp.savefig(transfer_experiments[i]['experiment_id'] + ".png", dpi=600, bbox_inches='tight')


def make_figure_5B():
    
    n_transfers = 1
    transfer_times = 30
    
    # parameters for figure 4C_1
    parameters = model_parameters.copy()
    parameters['transfer_times'] = n_transfers*[transfer_times]
    parameters['d_WT'] = 0
    parameters['d_LP'] = 0
    parameters['d_SP'] = 0
    
    # initial conditions for first transfer
    initial_host = model_parameters['initial_host']
    initial_phage = model_parameters['initial_MOI']*initial_host
    initial_condition_1 = [initial_host, 0, initial_phage, 0, 0]
    initial_condition_2 = [initial_host, 0, 0, initial_phage, 0]
    
    # this defines the specific transfer experiments we want to solve
    # each experiment consists of a set of biological/technical parameters (adsorption rate, burst size, transfer protocol, etc.),
    # and the initial condition for the very first transfer period
    transfer_experiments = list([
                                {
                                    'experiment_id': 'Fig_5B_1',
                                    'parameters': parameters,
                                    'initial_condition': initial_condition_1
                                },
                                {
                                    'experiment_id': 'Fig_5B_2',
                                    'parameters': parameters,
                                    'initial_condition': initial_condition_2
                                }
                            ])
    
    # run all defined transfer experiments in parallel
    solutions = run_parallel_experiments(transfer_experiments)

    
    fig, ax = pp.subplots(1, 1, figsize=(6,4))
    
    ax.set_xlabel('Time (minutes)')
    ax.set_ylabel('Phage density')
    ax.set_yscale('log')
    ax.set_ylim(5*10**3, 10**10)

    for i, solution in enumerate(solutions):
        sol_t = solution[:, 0]
        sol_WT = solution[:, 3]
        sol_LP = solution[:, 4]
        sol_SP = solution[:, 5]
        
        transfer_i = np.where(sol_t == 0)[0]
        transfer_i = np.append(transfer_i, len(sol_t)-1)
        
        plot_t = sol_t[transfer_i[0]:transfer_i[1]]
        
        ax.plot(plot_t, sol_WT[transfer_i[0]:transfer_i[1]], '-', lw=3, color='#29485d', label='ancestor (fast adsorption)')
        ax.plot(plot_t, sol_LP[transfer_i[0]:transfer_i[1]], '-', lw=3, color='#ae017e', label='slow adsorption')
        ax.plot(plot_t, sol_SP[transfer_i[0]:transfer_i[1]], '-', lw=3, color='crimson', label='fast adsorption')
  
    arrow_width = 12
    
    ax.annotate('', (16.9, 4*10**9), (0, 4*10**9), arrowprops=dict(arrowstyle='-',  capstyle='butt', facecolor='#bbd7eb', ec='#bbd7eb', linewidth=arrow_width))
    ax.text(4, 3.5*10**9, '1. adsorption phase', fontsize='x-small', fontweight='bold', color='dimgray')
  
    ax.annotate('', (21, 4*10**9), (16.8, 4*10**9), arrowprops=dict(arrowstyle='-',  capstyle='butt', facecolor='darkorange', ec='darkorange', linewidth=arrow_width))
    ax.text(18, 3.5*10**9, 'burst', fontsize='x-small', fontweight='bold', color='dimgray')
    
    ax.annotate('', (30, 4*10**9), (20.9, 4*10**9), arrowprops=dict(arrowstyle='-',  capstyle='butt', facecolor='#bbd7eb', ec='#bbd7eb', linewidth=arrow_width))
    ax.text(21.75, 3.5*10**9, '2. adsorption phase', fontsize='x-small', fontweight='bold', color='dimgray')
    
    ax.text(5, 10**7, 'slow adsorber', fontsize='small', fontweight='bold', color='#ae017e')
    ax.text(2, 10**5, 'fast adsorber', fontsize='small', fontweight='bold', color='#29485d')
  
    handles, labels = ax.get_legend_handles_labels()
    
    pp.tight_layout()
    pp.savefig("Fig_5B.png", dpi=600, bbox_inches='tight')


make_figure_4C()
make_figure_5B()
