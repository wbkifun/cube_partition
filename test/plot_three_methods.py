'''
Visualize the three-way comparison: SFC vs Stripe vs METIS
Generates publication-quality figures.
'''

from os.path import dirname, abspath, join
import sys
import time
import numpy as np
import pymetis

current_dir = dirname(abspath(__file__))
lib_dir = dirname(current_dir)
sys.path.append(lib_dir)

from cube_partition_sfc import CubePartitionSFC
from cube_partition_stripe import CubePartitionStripe
from cube_neighbor import CubeNeighbor
from compare_three_methods import build_cubed_sphere_adjacency, metis_partition, compute_metrics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def collect_data(ne, nproc_list, ngq=4):
    '''Collect comparison data for all nproc values.'''
    adjacency = build_cubed_sphere_adjacency(ne)

    sfc_pr, band_pr, metis_pr = [], [], []
    sfc_cr, band_cr, metis_cr = [], [], []
    sfc_tc, band_tc, metis_tc = [], [], []
    sfc_pr_std, band_pr_std, metis_pr_std = [], [], []

    for nproc in nproc_list:
        band_obj = CubePartitionStripe(ne, nproc)

        sfc_obj = CubePartitionSFC(ne, nproc)
        _, cr_sfc, _ = sfc_obj.make_cube_rank()
        m = compute_metrics(ne, nproc, cr_sfc, band_obj, ngq)
        sfc_pr.append(m['perimeter_ratio_mean'])
        sfc_pr_std.append(m['perimeter_ratio_std'])
        sfc_cr.append(m['comm_ratio_mean'])
        sfc_tc.append(m['total_comm'])

        _, cr_band, _ = band_obj.make_cube_rank()
        m = compute_metrics(ne, nproc, cr_band, band_obj, ngq)
        band_pr.append(m['perimeter_ratio_mean'])
        band_pr_std.append(m['perimeter_ratio_std'])
        band_cr.append(m['comm_ratio_mean'])
        band_tc.append(m['total_comm'])

        cr_metis = metis_partition(ne, nproc, adjacency)
        m = compute_metrics(ne, nproc, cr_metis, band_obj, ngq)
        metis_pr.append(m['perimeter_ratio_mean'])
        metis_pr_std.append(m['perimeter_ratio_std'])
        metis_cr.append(m['comm_ratio_mean'])
        metis_tc.append(m['total_comm'])

        print(f'  Ne={ne}, Nproc={nproc}: done')

    return {
        'nproc': np.array(nproc_list),
        'sfc': {'pr': np.array(sfc_pr), 'pr_std': np.array(sfc_pr_std),
                'cr': np.array(sfc_cr), 'tc': np.array(sfc_tc)},
        'stripe': {'pr': np.array(band_pr), 'pr_std': np.array(band_pr_std),
                   'cr': np.array(band_cr), 'tc': np.array(band_tc)},
        'metis': {'pr': np.array(metis_pr), 'pr_std': np.array(metis_pr_std),
                  'cr': np.array(metis_cr), 'tc': np.array(metis_tc)},
    }


def plot_comparison(data_list, ne_list, save=True):
    '''
    Plot three-panel figure:
      (a)(b)(c) Perimeter ratio for Ne=30, 60, 120
    '''
    fig, axes = plt.subplots(len(ne_list), 1, figsize=(8, 4 * len(ne_list)), sharex=False)
    if len(ne_list) == 1:
        axes = [axes]

    labels = ['(a)', '(b)', '(c)', '(d)']
    ylims = {30: 1.2, 60: 0.6, 120: 0.3}

    for idx, (data, ne) in enumerate(zip(data_list, ne_list)):
        ax = axes[idx]
        x = data['nproc']

        ax.plot(x, data['sfc']['pr'], 'bs-', ms=4, lw=1.2, label='SFC Partitioning')
        ax.fill_between(x,
                        data['sfc']['pr'] - data['sfc']['pr_std'],
                        data['sfc']['pr'] + data['sfc']['pr_std'],
                        alpha=0.15, color='b')

        ax.plot(x, data['metis']['pr'], 'g^-', ms=4, lw=1.2, label='METIS')
        ax.fill_between(x,
                        data['metis']['pr'] - data['metis']['pr_std'],
                        data['metis']['pr'] + data['metis']['pr_std'],
                        alpha=0.15, color='g')

        ax.plot(x, data['stripe']['pr'], 'ro-', ms=4, lw=1.2, label='Stripe Partitioning')
        ax.fill_between(x,
                        data['stripe']['pr'] - data['stripe']['pr_std'],
                        data['stripe']['pr'] + data['stripe']['pr_std'],
                        alpha=0.15, color='r')

        ax.set_title(f'{labels[idx]} $N_e$ = {ne}', fontsize=13)
        ax.set_ylabel('Mean Perimeter Ratio', fontsize=11)
        ax.set_ylim(0, ylims.get(ne, 1.0))
        ax.legend(loc='lower right', fontsize=9)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel('Number of Processes ($N_{proc}$)', fontsize=11)

    plt.tight_layout()
    if save:
        fname = f'png/three_methods_perimeter_ratio.png'
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        print(f'Saved: {fname}')
    plt.close()


def plot_traffic_reduction(data_list, ne_list, save=True):
    '''
    Plot total communication traffic for SFC, METIS, and Stripe.
    '''
    fig, ax = plt.subplots(1, 1, figsize=(9, 5.5))

    # ne_list is assumed [30, 60, 120] → linestyles: dotted, dashed, solid
    linestyles = {30: ':', 60: '--', 120: '-'}
    colors = {'sfc': 'b', 'metis': 'g', 'stripe': 'r'}
    markers = {'sfc': 's', 'metis': '^', 'stripe': 'o'}

    # Plot order: Stripe first, then METIS, then SFC (top to bottom in chart)
    # But build handles for custom legend order
    handles_map = {}

    for idx, (data, ne) in enumerate(zip(data_list, ne_list)):
        x = data['nproc']
        ls = linestyles[ne]

        for method, label_base in [('stripe', 'Stripe'), ('metis', 'METIS'), ('sfc', 'SFC')]:
            h, = ax.plot(x, data[method]['tc'], marker=markers[method],
                         linestyle=ls, color=colors[method],
                         ms=4, lw=1.2, label=f'{label_base} ($N_e$={ne})')
            handles_map[(method, ne)] = h

    ax.set_xlabel('Number of Processes ($N_{proc}$)', fontsize=11)
    ax.set_ylabel('Total Communication Traffic (points)', fontsize=11)
    ax.set_title('Spectral Element Method ($N_p=4$)', fontsize=13)
    ax.grid(True, alpha=0.3)

    # Custom legend: 3x3, rows=method, cols=Ne descending
    legend_order = []
    for method in ['stripe', 'metis', 'sfc']:
        for ne in reversed(ne_list):  # 120, 60, 30
            legend_order.append((method, ne))

    ordered_handles = [handles_map[k] for k in legend_order]
    ordered_labels = [h.get_label() for h in ordered_handles]
    ax.legend(ordered_handles, ordered_labels, loc='upper left', fontsize=8, ncol=3)

    plt.tight_layout()
    if save:
        fname = f'png/three_methods_traffic_reduction.png'
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        print(f'Saved: {fname}')
    plt.close()


def plot_summary_table(data_list, ne_list, save=True):
    '''
    Plot a summary bar chart comparing the three methods.
    '''
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ne_labels = [f'$N_e$={ne}' for ne in ne_list]
    x_pos = np.arange(len(ne_list))
    width = 0.25

    avg_pr_sfc = [d['sfc']['pr'].mean() for d in data_list]
    avg_pr_metis = [d['metis']['pr'].mean() for d in data_list]
    avg_pr_stripe = [d['stripe']['pr'].mean() for d in data_list]

    axes[0].bar(x_pos - width, avg_pr_sfc, width, label='SFC', color='steelblue')
    axes[0].bar(x_pos, avg_pr_metis, width, label='METIS', color='seagreen')
    axes[0].bar(x_pos + width, avg_pr_stripe, width, label='Stripe', color='indianred')
    axes[0].set_ylabel('Mean Perimeter Ratio (avg over all $N_{proc}$)')
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(ne_labels)
    axes[0].legend()
    axes[0].set_title('(a) Average Perimeter Ratio')
    axes[0].grid(True, alpha=0.3, axis='y')

    avg_tc_sfc = [d['sfc']['tc'].mean() for d in data_list]
    avg_tc_metis = [d['metis']['tc'].mean() for d in data_list]
    avg_tc_stripe = [d['stripe']['tc'].mean() for d in data_list]

    axes[1].bar(x_pos - width, avg_tc_sfc, width, label='SFC', color='steelblue')
    axes[1].bar(x_pos, avg_tc_metis, width, label='METIS', color='seagreen')
    axes[1].bar(x_pos + width, avg_tc_stripe, width, label='Stripe', color='indianred')
    axes[1].set_ylabel('Avg Total Communication Traffic')
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(ne_labels)
    axes[1].legend()
    axes[1].set_title('(b) Average Communication Traffic')
    axes[1].grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save:
        fname = f'png/three_methods_summary.png'
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        print(f'Saved: {fname}')
    plt.close()


if __name__ == '__main__':
    import os
    os.makedirs('png', exist_ok=True)

    ne_list = [30, 60, 120]
    data_list = []

    for ne in ne_list:
        total_elems = 6 * ne * ne
        max_nproc = min(total_elems // 4, 200)
        step = max(1, max_nproc // 15)
        nproc_list = [n for n in range(step, max_nproc + 1, step) if n >= 4]

        print(f'\nCollecting data for Ne={ne} ({len(nproc_list)} points)...')
        data = collect_data(ne, nproc_list)
        data_list.append(data)

    print('\nGenerating plots...')
    plot_comparison(data_list, ne_list)
    plot_traffic_reduction(data_list, ne_list)
    plot_summary_table(data_list, ne_list)
    print('Done!')
