"""Module with functions to carry out analyses for the TL project"""
from __future__ import division, with_statement
import partitions as parts
import numpy as np
from scipy import stats
import random
import csv
import signal
from contextlib import contextmanager

class TimeoutException(Exception): pass

@contextmanager
def time_limit(seconds):
    """Function to skip step after given time"""
    def signal_handler(signum, frame):
        raise TimeoutException, 'Time out!'
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

def get_QN_mean_var_data(data_dir):
    """Read in data file with study, Q, and N"""
    data = np.genfromtxt(data_dir, dtype = 'S25, i15, i15, f15, f15', delimiter = '\t', 
                         names = ['study', 'Q', 'N', 'mean', 'var'])
    return data

def get_study_info(data_dir):
    """Read in data file with study, taxon, and type"""
    data = np.genfromtxt(data_dir, dtype = 'S25, S25, S25', delimiter = '\t',
                          names = ['study', 'taxon', 'type'])
    return data

def RandomComposition_weak_xx(q, n):
    indices = sorted(np.random.randint(0, q, n - 1))
    parts = [(indices + [q])[i] - ([0] + indices)[i] for i in range(len(indices)+1)]
    return parts

def rand_compositions(q, n, sample_size, zeros):
    comps = []
    while len(comps) < sample_size:
        
        comp = RandomComposition_weak(q, n)
                
        if len(comp) != n or sum(comp) != q:
            print zeros,'error: q=',q,'!=sum(comp)=',sum(comp),'or n=',n,'!=len(comp)=',len(comp)
            sys.exit()
        comp.sort()
        comp.reverse()
        comps.append(comp)
    
    comps = [list(x) for x in set(tuple(x) for x in comps)]
    return comps

def get_var_for_Q_N(q, n, sample_size, t_limit, analysis):
    """Given q and n, returns a list of variance of length sample size with variance of 
    
    each sample partitions or compositions.
    
    """
    QN_var = []
    try:
        with time_limit(t_limit):
            for Niter in range(sample_size):
                if analysis == 'partition':
                    QN_parts = parts.rand_partitions(q, n, 1, 'bottom_up', {}, True)
                else: QN_parts = rand_compositions(q, n, 1, True)
                QN_var.append(np.var(QN_parts[0], ddof = 1))
            return QN_var
    except TimeoutException, msg:
        print 'Timed out!'
        return QN_var

def sample_var(data, study, sample_size = 1000, t_limit = 7200, analysis = 'partition'):
    """Obtain and record the variance of partition or composition samples.
    
    Input:
    data - data list read in with get_QN_mean_var_data()
    study - ID of study
    sample_size - number of samples to be drawn, default value is 1000
    t_limit - abort sampling procedure for one Q-N combo after t_limit seconds, default value is 7200 (2 hours)
    analysis - partition or composition
    
    """
    data_study = data[data['study'] == study]
    var_parts = []
    for record in data_study:
        q = record[1]
        n = record[2]
        out_row = [x for x in record]
        QN_var = get_var_for_Q_N(q, n, sample_size, t_limit, analysis)
        if len(QN_var) == sample_size:
            out_row.extend(QN_var)
            var_parts.append(out_row)
        else: break # Break out of for-loop if a Q-N combo is skipped
    
    if len(data_study) == len(var_parts): # If no QN combos are omitted, print to file
        out_write_var = open('taylor_QN_var_predicted_' + analysis + '_full.txt', 'a')
        for var_row in var_parts:
            print>>out_write_var, '\t'.join([str(x) for x in var_row])
        out_write_var.close()

def get_z_score(emp_var, sim_var_list):
    """Return the z-score as a measure of the discrepancy between empirical and sample variance"""
    sd_sim = (np.var(sim_var_list, ddof = 1)) ** 0.5
    return (emp_var - np.mean(sim_var_list)) / sd_sim
    
def TL_analysis(data, study, sample_size = 1000, t_limit = 7200, analysis = 'partition'):
    """Compare empirical TL relationship of one dataset to that obtained from random partitions or compositions."""
    data_study = data[data['study'] == study]
    data_study = data_study[data_study['N'] > 2] # Remove Q-N combos with N = 2
    var_parts = []
    for combo in data_study:
        q = combo[1]
        n = combo[2]
        QN_var = get_var_for_Q_N(q, n, sample_size, t_limit, analysis)
        if len(QN_var) == sample_size:
            var_parts.append(QN_var)
        else: break # Break out of for-loop if a Q-N combo is skipped
    
    if len(data_study) == len(var_parts): # IF no QN combos are omitted
        # 1. Predicted var for each Q-N combo
        var_study = np.zeros((len(data_study), ), dtype = [('f0', 'S25'), ('f1', int), ('f2', int), ('f3', float), 
                                                           ('f4', float), ('f5', float), ('f6', float)])
        var_study['f0'] = np.array([study] * len(data_study))
        var_study['f1'] = data_study['Q']
        var_study['f2'] = data_study['N']
        var_study['f3'] = np.array([np.mean(QN_var) for QN_var in var_parts])
        var_study['f4'] = np.array([np.median(QN_var) for QN_var in var_parts])
        var_study['f5'] = np.array([np.percentile(QN_var, 2.5) for QN_var in var_parts])
        var_study['f6'] = np.array([np.percentile(QN_var, 97.5) for QN_var in var_parts])
        out_write_var = open('taylor_QN_var_predicted_' + analysis + '.txt', 'a')
        out_var = csv.writer(out_write_var, delimiter = '\t')
        out_var.writerows(var_study)
        out_write_var.close()
        
        # 2. Predicted form of TL for study
        b_list = []
        inter_list = []
        psig = 0
        R2_list = []
        effective_sample = 0
        for i in range(sample_size):
            var_list = np.array([var_part[i] for var_part in var_parts])
            mean_list = data_study['mean']
            mean_list = mean_list[var_list != 0] # Omit samples of zero variance in computing TL
            var_list = var_list[var_list != 0]
            b, inter, rval, pval, std_err = stats.linregress(np.log(mean_list), np.log(var_list))
            b_list.append(b)
            inter_list.append(inter)
            R2_list.append(rval ** 2)
            if pval < 0.05: psig += 1
        psig = psig / sample_size
        OUT = open('taylor_form_predicted_' + analysis + '.txt', 'a')
        print>>OUT, study, psig, np.mean(R2_list), np.mean(b_list), np.median(b), np.percentile(b_list, 2.5), \
             np.percentile(b_list, 97.5), np.mean(inter_list), np.median(inter_list), np.percentile(inter_list, 2.5), \
             np.percentile(inter_list, 97.5)
        OUT.close()

def inclusion_criteria(dat_study, sig = False):
    """Criteria that datasets need to meet to be included in the analysis"""
    b, inter, rval, pval, std_err = stats.linregress(np.log(dat_study['mean']), np.log(dat_study['var']))
    dat_study = dat_study[dat_study['N'] > 2] # Doesn't make too much sense to talk about variance among two points
    if len(dat_study) >= 5: # More than 5 observations
        if len(dat_study[dat_study['Q'] > 5]) / len(dat_study) >= 0.5: # Not predominantly small values
            if ((not sig) or (pval < 0.05)): # If significance is not required, or if the relationship is significant
                return True
    else: return False