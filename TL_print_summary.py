"""Print the summary for Results"""
from __future__ import division
import TL_functions as tl
import numpy as np
from scipy import stats

study_info = tl.get_study_info('study_taxon_type.txt')
tl_pars_par = tl.get_tl_par_file('TL_form_partition.txt')
tl_pars_comp = tl.get_tl_par_file('TL_form_composition.txt')

var_par = tl.get_var_sample_file('taylor_QN_var_predicted_partition_1000_full.txt')
var_comp = tl.get_var_sample_file('taylor_QN_var_predicted_composition_1000_full.txt')

study_spatial = [study for study in np.unique(var_par['study']) if study_info['type'] [study_info['study'] == study]== 'spatial']
study_temporal = [study for study in np.unique(var_par['study']) if study_info['type'][study_info['study'] == study] == 'temporal']
# 1. Curvature
par_quad = tl.get_val_ind_sample_file('TL_quad_p_partition.txt')
comp_quad = tl.get_val_ind_sample_file('TL_quad_p_composition.txt')

sig_spatial, sig_temporal, sig_par, sig_comp, tot_sig = 0, 0, 0, 0, 0
for study in study_spatial:
    row_study_par = list(par_quad[par_quad['study'] == study][0])
    row_study_comp = list(comp_quad[comp_quad['study'] == study][0])
    if row_study_par[1] < 0.05: sig_spatial += 1
    sig_par += len([x for x in row_study_par[2:] if x < 0.05])
    sig_comp += len([x for x in row_study_comp[2:] if x < 0.05])
    tot_sig += len(row_study_comp[2:])
    
for study in study_temporal:
    row_study_par = list(par_quad[par_quad['study'] == study][0])
    row_study_comp = list(comp_quad[comp_quad['study'] == study][0])
    if row_study_par[1] < 0.05: sig_temporal += 1
    sig_par += len([x for x in row_study_par[2:] if x < 0.05])
    sig_comp += len([x for x in row_study_comp[2:] if x < 0.05])
    tot_sig += len(row_study_comp[2:])

print "Number of spatial TLs with curvature: ", str(sig_spatial)
print "Number of temporal TLs with curvature: ", str(sig_temporal)
print "Proportion of partitions with curvature: ", str(sig_par / tot_sig)
print "Proportion of compositions with curvature: ", str(sig_comp / tot_sig)

# 2. Exponent b, p, avg r^2
b_obs, b_par, b_comp, b_type = [], [], [], []
r2_obs, r2_par, r2_comp, p_obs, p_par, p_comp = [], [], [], [], [], []
for study in np.unique(var_par['study']):
    b_obs.append((tl_pars_par['b_obs'][tl_pars_par['study'] == study])[0])
    b_type.append((study_info['type'][study_info['study'] == study])[0])
    r2_obs.append((tl_pars_par['R2_obs'][tl_pars_par['study'] == study])[0])
    p_obs.append((tl_pars_par['p_obs'][tl_pars_par['study'] == study])[0])
    sample_par = var_par[var_par['study'] == study]
    sample_comp = var_comp[var_comp['study'] == study]
    for i in xrange(1000):
        sample_par_i = [sample_par[x][i + 5] for x in xrange(len(sample_par))]
        mean_par = [sample_par['mean'][p] for p in xrange(len(sample_par)) if sample_par_i[p] != 0]
        sample_par_i = [sample_par_i[p] for p in xrange(len(sample_par_i)) if sample_par_i[p] != 0]
        b, inter, r, p, std_error = stats.linregress(np.log(mean_par), np.log(sample_par_i))
        b_par.append(b)
        r2_par.append(r ** 2)
        p_par.append(p)
        sample_comp_i = [sample_comp[y][i + 5] for y in xrange(len(sample_comp))]
        mean_comp = [sample_comp['mean'][q] for q in xrange(len(sample_comp)) if sample_comp_i[q] != 0]
        sample_comp_i = [sample_comp_i[q] for q in xrange(len(sample_comp_i)) if sample_comp_i[q] != 0]        
        b, inter, r, p, std_error = stats.linregress(np.log(mean_comp), np.log(sample_comp_i))
        b_comp.append(b)
        r2_comp.append(r ** 2)
        p_comp.append(p)

b_spatial = [b_obs[i] for i in range(len(b_obs)) if b_type[i] == 'spatial']
b_temporal = [b_obs[i] for i in range(len(b_obs)) if b_type[i] == 'temporal']
r2_spatial = [r2_obs[i] for i in range(len(r2_obs)) if b_type[i] == 'spatial']
r2_temporal = [r2_obs[i] for i in range(len(r2_obs)) if b_type[i] == 'temporal']
p_spatial = [p_obs[i] for i in range(len(p_obs)) if b_type[i] == 'spatial']
p_temporal = [p_obs[i] for i in range(len(p_obs)) if b_type[i] == 'temporal']
print "b in partitions - min, max, proportion between 1 & 2: ", str(min(b_par)), " ,", \
      str(max(b_par)), " ,", str(len([x for x in b_par if 1 < x < 2]) / len(b_par))
print "b in compositions - min, max, proportion between 1 & 2: ", str(min(b_comp)), " ,", \
      str(max(b_comp)), " ,", str(len([x for x in b_comp if 1 < x < 2]) / len(b_comp))
print "b in spatial TL - min, max, proportion between 1 & 2: ", str(min(b_spatial)), " ,", \
      str(max(b_spatial)), " ,", str(len([x for x in b_spatial if 1 < x < 2]) / len(b_spatial))
print "b in temporal TL - min, max, proportion between 1 & 2: ", str(min(b_temporal)), " ,", \
      str(max(b_temporal)), " ,", str(len([x for x in b_temporal if 1 < x < 2]) / len(b_temporal))
print "Proportion of TL from partitions that are significant: ", str(len([p for p in p_par if p < 0.05]) / 1000 / 111)
print "Proportion of TL from compositions that are significant: ", str(len([p for p in p_comp if p < 0.05]) / 1000 / 111)
print "Proportion of spatial TL that are significant: ", str(len([p for p in p_spatial if p < 0.05]) / len(p_spatial))
print "Proportion of temporal TL that are significant: ", str(len([p for p in p_temporal if p < 0.05]) / len(p_temporal))
print "R2 from partitions - min, max, average: ", str(min(r2_par)), " ,", str(max(r2_par)), " ,", str(np.mean(r2_par))
print "R2 from compositios - min, max, average: ", str(min(r2_comp)), " ,", str(max(r2_comp)), " ,", str(np.mean(r2_comp))
print "R2 from spatial TL - min, max, average: ", str(min(r2_spatial)), " ,", str(max(r2_spatial)), " ,", str(np.mean(r2_spatial))
print "R2 from temporal TL - min, max, average: ", str(min(r2_temporal)), " ,", str(max(r2_temporal)), " ,", str(np.mean(r2_temporal))

# 3. var and b from the feasible set against emp values
var_out_spa_par, var_out_spa_comp, var_out_temp_par, var_out_temp_comp = 0, 0, 0, 0
var_tot_spa, var_tot_temp = 0, 0
b_out_spa_par, b_out_spa_comp, b_out_temp_par, b_out_temp_comp = 0, 0, 0, 0

for study in study_sig_spatial:
    var_par_study = var_par[var_par['study'] == study]
    var_tot_spa += len(var_par_study)
    for row in var_par_study: 
        var_feas = list(row)[5:]
        if not np.percentile(var_feas, 2.5) < row[4] < np.percentile(var_feas, 97.5): var_out_spa_par += 1
    var_comp_study = var_comp[var_comp['study'] == study]
    for row in var_comp_study: 
        var_feas = list(row)[5:]
        if not np.percentile(var_feas, 2.5) < row[4] < np.percentile(var_feas, 97.5): var_out_spa_comp += 1
    
    b_row_par = tl_pars_par[tl_pars_par['study'] == study]
    if not b_row_par['b_lower'] < b_row_par['b_obs'] < b_row_par['b_upper']: b_out_spa_par += 1
    b_row_comp = tl_pars_comp[tl_pars_comp['study'] == study]
    if not b_row_comp['b_lower'] < b_row_comp['b_obs'] < b_row_comp['b_upper']: b_out_spa_comp += 1
    
for study in study_sig_temporal:
    var_par_study = var_par[var_par['study'] == study]
    var_tot_temp += len(var_par_study)
    for row in var_par_study: 
        var_feas = list(row)[5:]
        if not np.percentile(var_feas, 2.5) < row[4] < np.percentile(var_feas, 97.5): var_out_temp_par += 1
    var_comp_study = var_comp[var_comp['study'] == study]
    for row in var_comp_study: 
        var_feas = list(row)[5:]
        if not np.percentile(var_feas, 2.5) < row[4] < np.percentile(var_feas, 97.5): var_out_temp_comp += 1
    
    b_row_par = tl_pars_par[tl_pars_par['study'] == study]
    if not b_row_par['b_lower'] < b_row_par['b_obs'] < b_row_par['b_upper']: b_out_temp_par += 1
    b_row_comp = tl_pars_comp[tl_pars_comp['study'] == study]
    if not b_row_comp['b_lower'] < b_row_comp['b_obs'] < b_row_comp['b_upper']: b_out_temp_comp += 1

print "Proportion of variance out of 95% quantile: spaital (par, comp), temporal (par, comp): ", \
      str(var_out_spa_par / var_tot_spa), str(var_out_spa_comp / var_tot_spa), str(var_out_temp_par / var_tot_temp), str(var_out_temp_comp / var_tot_temp)
print "Proortion of variance out of 95% quantile: partition, composition: ", \
      str((var_out_spa_par + var_out_temp_par) / (var_tot_spa + var_tot_temp)), str((var_out_spa_comp + var_out_temp_comp) / (var_tot_spa + var_tot_temp))
print "Number of b out of 95% quantile: spatial (par, comp), temporal (par, comp):", \
      str(b_out_spa_par), str(b_out_spa_comp), str(b_out_temp_par), str(b_out_temp_comp)
