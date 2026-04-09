import numpy as np
from matplotlib import pyplot as pp
from matplotlib import colormaps as cmaps
from model import model_parameters, run_parallel_experiments


def make_figure_4C():

    n_transfers = 20
    transfer_times = 180

    # parameters for figure
    parameters_adsorption_only = model_parameters.copy()
    parameters_adsorption_only['transfer_times'] = n_transfers*[transfer_times]
    parameters_adsorption_only['d_WT'] = 0
    parameters_adsorption_only['d_LP'] = 0
    parameters_adsorption_only['d_SP'] = 0

    # parameters for figure
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

    # parameters for figure
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


def make_figure_6BC():

    n_transfers = 20
    transfer_times = 30

    # parameters for figure
    parameters = model_parameters.copy()
    parameters['transfer_times'] = n_transfers*[transfer_times]

    # initial conditions for first transfer
    initial_host = model_parameters['initial_host']
    initial_phage = model_parameters['initial_MOI']*initial_host
    initial_condition = [initial_host, 0, initial_phage, 0, 0]

    transfer_experiments = list([
                                {
                                    'experiment_id': 'Fig_6BC',
                                    'parameters': parameters,
                                    'initial_condition': initial_condition
                                }
                            ])

    # run all defined transfer experiments in parallel
    solutions = run_parallel_experiments(transfer_experiments)


    for i, solution in enumerate(solutions):
        sol_t = solution[:, 0]
        sol_B = solution[:, 1]
        sol_Bi = solution[:, 2]
        sol_WT = solution[:, 3]

        transfer_i = np.where(sol_t == 0)[0]
        transfer_i = np.append(transfer_i, len(sol_t)-1)

        fig, axs = pp.subplots(1, 2, sharex=False, figsize=(8, 3.5))

        axs[0].set_xlabel('Time (mins)')
        axs[0].set_ylabel('Phage density')
        axs[0].set_yscale('log')
        axs[0].set_ylim(5*10**3, 10**12)

        axs[1].set_xlabel('Time (mins)')
        axs[1].set_ylabel('Bacterial density')
        axs[1].set_ylim(10**5, 5*10**8)
        axs[1].set_yscale('log')

        n_plot_seasons = len(transfer_i)-1

        greys = cmaps['Greys'](np.linspace(0.05, 1, n_plot_seasons))
        blues = cmaps['Blues'](np.linspace(0.05, 1, n_plot_seasons))
        oranges = cmaps['Oranges'](np.linspace(0.05, 1, n_plot_seasons))

        for j in range(10, n_plot_seasons):

            axs[0].plot(sol_t[transfer_i[j]:transfer_i[j+1]], sol_WT[transfer_i[j]:transfer_i[j+1]], '-', lw=2, color=greys[j])

            axs[1].plot(sol_t[transfer_i[j]:transfer_i[j+1]], sol_B[transfer_i[j]:transfer_i[j+1]], '-', lw=2, color=blues[j], label='not infected')
            axs[1].plot(sol_t[transfer_i[j]:transfer_i[j+1]], sol_Bi[transfer_i[j]:transfer_i[j+1]], '-', lw=2, color=oranges[j], label='infected')


        handles, labels = axs[1].get_legend_handles_labels()
        axs[1].legend(handles=handles[-2:], labels=labels[-2:], ncols=2, frameon=False)

        arrow_width = 8

        axs[0].annotate('', (17, 4*10**11), (0, 4*10**11), arrowprops=dict(facecolor='#bbd7eb', ec='#bbd7eb', width=arrow_width, headwidth=1.5*arrow_width))
        axs[0].text(4, 3.4*10**11, '1st phage cycle', fontsize='xx-small', fontweight='bold', color=greys[int(0.75*len(greys))])

        axs[0].annotate('', (30.3, 10**11), (16.7, 10**11), arrowprops=dict(arrowstyle='-',  capstyle='butt', facecolor='#bbd7eb', ec='#bbd7eb', linewidth=arrow_width))
        axs[0].text(19, 8.5*10**10, '2nd phage cycle', fontsize='xx-small', fontweight='bold', color=greys[int(0.75*len(greys))])

        axs[0].annotate('', (4, 10**11), (0, 10**11), arrowprops=dict(facecolor='#bbd7eb', ec='#bbd7eb', width=arrow_width))

        axs[0].annotate('', (23, 2.9*10**10), (4, 2.9*10**10), arrowprops=dict(facecolor='#bbd7eb', ec='#bbd7eb', width=arrow_width, headwidth=1.5*arrow_width))
        axs[0].text(9, 2.5*10**10, '3rd phage cycle', fontsize='xx-small', fontweight='bold', color=greys[int(0.75*len(greys))])

        pp.tight_layout()
        pp.savefig("Fig_6BC.png", dpi=600, bbox_inches='tight')


make_figure_4C()
make_figure_5B()
make_figure_6BC()
