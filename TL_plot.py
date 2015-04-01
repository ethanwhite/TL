## Test the plotting function
from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import TL_functions as tl
import numpy as np
from scipy import stats
import random

# Figure 1 - visual representation using three studies
study_list = ['1_1', '10_1', '52_11']
fig = plt.figure(figsize = (10.5, 7))
iplot = 1
for feas_type in ['partition', 'composition']:
    for study in study_list:
        ax = plt.subplot(2, 3, iplot)
        if iplot == 1 or iplot == 4: legend = True
        else: legend = False
        tl.plot_emp_vs_sim(study, feas_type = feas_type, ax = ax, legend = legend)     
        iplot += 1
plt.subplots_adjust(wspace = 0.29, hspace = 0.29)
plt.savefig('Fig1.pdf', dpi = 600)

# Figure 2 - compare the full distribution of empirical TLs and those from the feasible sets
study_info = tl.get_study_info('study_taxon_type.txt')
tl_pars_par = tl.get_tl_par_file('out_files/TL_form_partition.txt')

var_par = tl.get_var_sample_file('out_files/taylor_QN_var_predicted_partition_1000_full.txt')
var_comp = tl.get_var_sample_file('out_files/taylor_QN_var_predicted_composition_1000_full.txt')
par_quad = tl.get_val_ind_sample_file('out_files/TL_quad_p_partition.txt')
comp_quad = tl.get_val_ind_sample_file('out_files/TL_quad_p_composition.txt')

b_obs, b_par, b_comp, b_type = [], [], [], []
p_obs, p_par, p_comp = [], [], []
pcurv_obs, pcurv_par, pcurv_comp = [], [], []
r2_obs, r2_par, r2_comp = [], [], []
for study in np.unique(var_par['study']):
    b_obs.append((tl_pars_par['b_obs'][tl_pars_par['study'] == study])[0])
    p_obs.append((tl_pars_par['p_obs'][tl_pars_par['study'] == study])[0])
    r2_obs.append((tl_pars_par['R2_obs'][tl_pars_par['study'] == study])[0])
    pcurv_obs.append((par_quad['emp_val'][par_quad['study'] == study])[0])
    b_type.append((study_info['type'][study_info['study'] == study])[0])
    
    sample_par = var_par[var_par['study'] == study]
    sample_comp = var_comp[var_comp['study'] == study]
    for i in xrange(1000):
        sample_par_i = [sample_par[x][i + 5] for x in xrange(len(sample_par))]
        mean_par = [sample_par['mean'][p] for p in xrange(len(sample_par)) if sample_par_i[p] != 0]
        sample_par_i = [sample_par_i[p] for p in xrange(len(sample_par_i)) if sample_par_i[p] != 0]
        b_par_i, inter_par_i, r_par_i, p_par_i, std_err_par_i = stats.linregress(np.log(mean_par), np.log(sample_par_i))
        b_par.append(b_par_i)
        p_par.append(p_par_i)
        r2_par.append(r_par_i ** 2)
        pcurv_par.append(par_quad[par_quad['study'] == study][0][i + 2])
        
        sample_comp_i = [sample_comp[y][i + 5] for y in xrange(len(sample_comp))]
        mean_comp = [sample_comp['mean'][q] for q in xrange(len(sample_comp)) if sample_comp_i[q] != 0]
        sample_comp_i = [sample_comp_i[q] for q in xrange(len(sample_comp_i)) if sample_comp_i[q] != 0] 
        b_comp_i, inter_comp_i, r_comp_i, p_comp_i, std_err_comp_i = stats.linregress(np.log(mean_comp), np.log(sample_comp_i))
        b_comp.append(b_comp_i)
        p_comp.append(p_comp_i)
        r2_comp.append(r_comp_i ** 2)
        pcurv_comp.append(comp_quad[comp_quad['study'] == study][0][i + 2])
        
fig = plt.figure(figsize = (7, 7)) 
ax_p = plt.subplot(221)
tl.plot_dens_par_comp(p_obs, p_par, p_comp, ax = ax_p, legend = True, loc = 1, vline = 0.05, xlim = [0, 0.2])
ax_p.annotate('(A)', xy = (0.05, 0.92), xycoords = 'axes fraction', fontsize = 10)
plt.xlabel('p-value for b', fontsize = 8)
plt.ylabel('Density', fontsize = 8)

ax_curv = plt.subplot(222)
tl.plot_dens_par_comp(pcurv_obs, pcurv_par, pcurv_comp, ax = ax_curv, vline = 0.05, xlim = [0, 1])
ax_curv.annotate('(B)', xy = (0.05, 0.92), xycoords = 'axes fraction', fontsize = 10)
plt.xlabel('p-value for quadratic term', fontsize = 8)
plt.ylabel('Density', fontsize = 8)

ax_r2 = plt.subplot(223)
tl.plot_dens_par_comp(r2_obs, r2_par, r2_comp, ax = ax_r2, xlim = [0, 1])
ax_r2.annotate('(C)', xy = (0.05, 0.92), xycoords = 'axes fraction', fontsize = 10)
plt.xlabel(r'$R^2$', fontsize = 8)
plt.ylabel('Density', fontsize = 8)

ax_b = plt.subplot(224) 
tl.plot_dens_par_comp(b_obs, b_par, b_comp, ax = ax_b, xlim = [0, 4])
ax_b.annotate('(D)', xy = (0.05, 0.92), xycoords = 'axes fraction', fontsize = 10)
plt.xlabel('Exponent b', fontsize = 8)
plt.ylabel('Density', fontsize = 8)
plt.savefig('Fig2.pdf', dpi = 600)

# Figure 3 - compare empirical versus feasible set variance and b
tl_pars_par = tl.get_tl_par_file('out_files/TL_form_partition.txt')
study_sig = tl_pars_par['study']

var_par = tl.get_var_sample_file('out_files/taylor_QN_var_predicted_partition_1000_full.txt')
var_comp = tl.get_var_sample_file('out_files/taylor_QN_var_predicted_composition_1000_full.txt')
par_index = [i for i in range(len(var_par)) if var_par['study'][i] in study_sig]
var_par = var_par[par_index]
comp_index = [i for i in range(len(var_comp)) if var_comp['study'][i] in study_sig]
var_comp = var_comp[comp_index]

study_info = tl.get_study_info('study_taxon_type.txt')
# Here the values are relative to the emp value
expc_upper_par, expc_upper_comp, expc_par, expc_comp, expc_lower_par, expc_lower_comp = [], [], [], [], [], []

for i in xrange(len(var_par)):
    expc_sample_par = [var_par[i][j] for j in xrange(5, 1005)]
    expc_par.append(np.mean(expc_sample_par))
    expc_lower_par.append(np.percentile(expc_sample_par, 2.5))
    expc_upper_par.append(np.percentile(expc_sample_par, 97.5))

    expc_sample_comp = [var_comp[i][j] for j in xrange(5, 1005)]
    expc_comp.append(np.mean(expc_sample_comp))
    expc_lower_comp.append(np.percentile(expc_sample_comp, 2.5) )
    expc_upper_comp.append(np.percentile(expc_sample_comp, 97.5))
    
fig = plt.figure(figsize = (7, 7))
ax_par = plt.subplot(221)
tl.plot_obs_expc_new(var_par['var'], expc_par, expc_upper_par, expc_lower_par, 'partition', True, ax = ax_par)
plt.xlabel(r'Index for  $s^2$', fontsize = 10)
plt.ylabel(r'$s_{partition}^2$ / $s_{empirical}^2$', fontsize = 12)
plt.title('Partitions')

ax_comp = plt.subplot(222)
tl.plot_obs_expc_new(var_comp['var'], expc_comp, expc_upper_comp, expc_lower_comp, 'composition', True, ax = ax_comp)
plt.xlabel(r'Index for  $s^2$', fontsize = 10)
plt.ylabel(r'$s_{composition}^2$/ $s_{empirical}^2$', fontsize = 12)
plt.title('Compositions')

ax_b_par = plt.subplot(223)
tl_pars_par = tl.get_tl_par_file('out_files/TL_form_partition.txt')
tl.plot_obs_expc_new(tl_pars_par['b_obs'], tl_pars_par['b_expc'], tl_pars_par['b_upper'], \
                     tl_pars_par['b_lower'], 'partition', False, ax = ax_b_par)
plt.xlabel('Index for b', fontsize = 10)
plt.ylabel(r'$b_{partition}$ - $b_{empirical}$', fontsize = 12)

ax_b_comp = plt.subplot(224)
tl_pars_comp = tl.get_tl_par_file('out_files/TL_form_composition.txt')
tl.plot_obs_expc_new(tl_pars_comp['b_obs'], tl_pars_comp['b_expc'], tl_pars_comp['b_upper'], \
                     tl_pars_comp['b_lower'], 'composition', False, ax = ax_b_comp)
plt.xlabel('Index for b', fontsize = 10)
plt.ylabel(r'$b_{composition}$ - $b_{empirical}$', fontsize = 12)

plt.subplots_adjust(wspace = 0.29, hspace = 0.29)
plt.savefig('Fig3.pdf', dpi = 600)

# 10. Figure B1 - examples of empirical variance versus the full distribution from the feasible set
var_par = tl.get_var_sample_file('out_files/taylor_QN_var_predicted_partition_1000_full.txt')
var_comp = tl.get_var_sample_file('out_files/taylor_QN_var_predicted_composition_1000_full.txt')
tl_pars_par = tl.get_tl_par_file('out_files/TL_form_partition.txt')
random.seed(4) 
qn_sets = random.sample(range(len(var_par)), 3) 

fig = plt.figure(figsize = (10.5, 3.5)) 
# Plot the 3 Q-N pairs 
for i, qn_pair in enumerate(qn_sets):
    ax_qn = plt.subplot(1, 3, i + 1)
    dat_row_par = list(var_par[qn_pair])
    study_comp = var_comp[var_comp['study'] == dat_row_par[0]]
    dat_row_comp = study_comp[(study_comp['Q'] == dat_row_par[1]) * (study_comp['N'] == dat_row_par[2])]
    dat_row_comp = list(dat_row_comp[0])
    if i == 0:
        tl.plot_dens_par_comp_single_obs(dat_row_par[4], dat_row_par[5:], dat_row_comp[5:], \
                                         ax = ax_qn, legend = True, loc = 2)
    else: tl.plot_dens_par_comp_single_obs(dat_row_par[4], dat_row_par[5:], dat_row_comp[5:], ax = ax_qn)
    ax_qn.annotate('Q = ' + str(dat_row_par[1]), xy = (0.7, 0.92), xycoords = 'axes fraction', fontsize = 10)
    ax_qn.annotate('N = ' + str(dat_row_par[2]), xy = (0.7, 0.82), xycoords = 'axes fraction', fontsize = 10)
    plt.xlabel('Variance', fontsize = 10)
    plt.ylabel('Density', fontsize = 10)
plt.subplots_adjust(wspace = 0.29)
plt.savefig('FigB1.pdf', dpi = 600)

# Figure B2 - results from 4000 samples
study_info = tl.get_study_info('study_taxon_type.txt')
tl_pars_par = tl.get_tl_par_file('out_files/TL_form_partition_4000.txt')

var_par_1000 = tl.get_var_sample_file('out_files/taylor_QN_var_predicted_partition_1000_full.txt')
var_par = tl.get_var_sample_file('out_files/taylor_QN_var_predicted_partition_4000_full.txt', sample_size = 4000)
var_comp = tl.get_var_sample_file('out_files/taylor_QN_var_predicted_composition_4000_full.txt', sample_size = 4000)
par_quad = tl.get_val_ind_sample_file('out_files/TL_quad_p_partition_4000.txt', sample_size = 4000)
comp_quad = tl.get_val_ind_sample_file('out_files/TL_quad_p_composition_4000.txt', sample_size = 4000)

b_obs, b_par, b_comp, b_type = [], [], [], []
p_obs, p_par, p_comp = [], [], []
pcurv_obs, pcurv_par, pcurv_comp = [], [], []
r2_obs, r2_par, r2_comp = [], [], []
for study in np.unique(var_par_1000['study']):
    b_obs.append((tl_pars_par['b_obs'][tl_pars_par['study'] == study])[0])
    p_obs.append((tl_pars_par['p_obs'][tl_pars_par['study'] == study])[0])
    r2_obs.append((tl_pars_par['R2_obs'][tl_pars_par['study'] == study])[0])
    pcurv_obs.append((par_quad['emp_val'][par_quad['study'] == study])[0])
    b_type.append((study_info['type'][study_info['study'] == study])[0])
    
    sample_par = var_par[var_par['study'] == study]
    sample_comp = var_comp[var_comp['study'] == study]
    for i in xrange(4000):
        sample_par_i = [sample_par[x][i + 5] for x in xrange(len(sample_par))]
        mean_par = [sample_par['mean'][p] for p in xrange(len(sample_par)) if sample_par_i[p] != 0]
        sample_par_i = [sample_par_i[p] for p in xrange(len(sample_par_i)) if sample_par_i[p] != 0]
        b_par_i, inter_par_i, r_par_i, p_par_i, std_err_par_i = stats.linregress(np.log(mean_par), np.log(sample_par_i))
        b_par.append(b_par_i)
        p_par.append(p_par_i)
        r2_par.append(r_par_i ** 2)
        pcurv_par.append(par_quad[par_quad['study'] == study][0][i + 2])
        
        sample_comp_i = [sample_comp[y][i + 5] for y in xrange(len(sample_comp))]
        mean_comp = [sample_comp['mean'][q] for q in xrange(len(sample_comp)) if sample_comp_i[q] != 0]
        sample_comp_i = [sample_comp_i[q] for q in xrange(len(sample_comp_i)) if sample_comp_i[q] != 0] 
        b_comp_i, inter_comp_i, r_comp_i, p_comp_i, std_err_comp_i = stats.linregress(np.log(mean_comp), np.log(sample_comp_i))
        b_comp.append(b_comp_i)
        p_comp.append(p_comp_i)
        r2_comp.append(r_comp_i ** 2)
        pcurv_comp.append(comp_quad[comp_quad['study'] == study][0][i + 2])
        
fig = plt.figure(figsize = (7, 7)) 
ax_p = plt.subplot(221)
tl.plot_dens_par_comp(p_obs, p_par, p_comp, ax = ax_p, legend = True, loc = 1, vline = 0.05, xlim = [0, 0.2])
ax_p.annotate('(A)', xy = (0.05, 0.92), xycoords = 'axes fraction', fontsize = 10)
plt.xlabel('p-value for b', fontsize = 8)
plt.ylabel('Density', fontsize = 8)

ax_curv = plt.subplot(222)
pcurv_par = [x for x in pcurv_par if not np.isnan(x)]
pcurv_comp = [x for x in pcurv_comp if not np.isnan(x)]
tl.plot_dens_par_comp(pcurv_obs, pcurv_par,  pcurv_comp, ax = ax_curv, vline = 0.05, xlim = [0, 1])
ax_curv.annotate('(B)', xy = (0.05, 0.92), xycoords = 'axes fraction', fontsize = 10)
plt.xlabel('p-value for quadratic term', fontsize = 8)
plt.ylabel('Density', fontsize = 8)

ax_r2 = plt.subplot(223)
tl.plot_dens_par_comp(r2_obs, r2_par, r2_comp, ax = ax_r2, xlim = [0, 1])
ax_r2.annotate('(C)', xy = (0.05, 0.92), xycoords = 'axes fraction', fontsize = 10)
plt.xlabel(r'$R^2$', fontsize = 8)
plt.ylabel('Density', fontsize = 8)

ax_b = plt.subplot(224) 
tl.plot_dens_par_comp(b_obs, b_par, b_comp, ax = ax_b, xlim = [0, 4])
ax_b.annotate('(D)', xy = (0.05, 0.92), xycoords = 'axes fraction', fontsize = 10)
plt.xlabel('Exponent b', fontsize = 8)
plt.ylabel('Density', fontsize = 8)
plt.savefig('FigB2.pdf', dpi = 600)
