from __future__ import division
import TL_functions as tl
import multiprocessing

def get_good_study(data_list, sig = True):
    """Return a list of studies that pass the criteria check"""
    study_list_data = sorted(list(set(data_list['study'])))
    good_study_list = []
    for study in study_list_data:
        data_study = data_list[data_list['study'] == study]
        if tl.inclusion_criteria(data_study, sig = sig):
            good_study_list.append(study)
    return good_study_list

data_lit = tl.get_QN_mean_var_data('data_literature.txt')
data_glenda = tl.get_QN_mean_var_data('data_Glenda.txt')
good_list_lit = get_good_study(data_lit, sig = False)
good_list_glenda = get_good_study(data_glenda)

def map_lit_partition(study):
    tl.TL_analysis(data_lit, study)

def map_lit_composition(study):
    tl.TL_analysis(data_lit, study, analysis = 'composition')

def map_glenda_partition(study):
    tl.TL_analysis(data_glenda, study)

def map_glenda_composition(study):
    tl.TL_analysis(data_glenda, study, analysis = 'composition')
      
step_list = [[map_lit_partition, good_list_lit], [map_glenda_partition, good_list_glenda], 
             [map_lit_composition, good_list_lit], [map_glenda_composition, good_list_glenda]]
for i in range(len(step_list)):
    step_to_take = step_list[i]
    pool = multiprocessing.Pool(8)
    pool.map(step_to_take[0], step_to_take[1])
    pool.close()
    pool.join()