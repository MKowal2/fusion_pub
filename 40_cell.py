import sys
import os
import csv
from os import listdir
from os.path import isfile, join
from glob import glob
import numpy as np
from itertools import combinations
import ast
from output_tool_grab import *
from evaluations import *

def create_40_cell_individual_cell_lines():
    root_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/labs'
    save_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/'

    lab_directories = [join(root_dir, f) for f in listdir(root_dir)]
    csv_fusion_master_list = [['dataset', '#_fusions_predicted', 'fusion_pred_list']]
    for lab_dir in lab_directories:
        lab = lab_dir.split('/')[-1]
        csv_list = [['dataset', 'cell_line', '#_fusions_predicted', 'fusion_pred_list']]
        master_fusion_list = []
        starfusion_pred = get_starfusion_dict(lab_dir+'*/', 'star-fusion.fusion_predictions.abridged.tsv')
        for cell in starfusion_pred[lab]:
            fusion_pred_list = starfusion_pred[lab][cell][0]
            new_list = [lab, cell, len(fusion_pred_list), fusion_pred_list]
            csv_list.append(new_list)

            master_fusion_list += fusion_pred_list

        new_master_list = [lab, len(master_fusion_list), master_fusion_list]
        csv_fusion_master_list.append(new_master_list)

        with open(save_dir + lab + '_starfusion_output'+ '.csv', "w", newline="") as f2:
            writer = csv.writer(f2)
            writer.writerows(csv_list)

    with open(save_dir + 'starfusion_individual_master.csv', "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(csv_fusion_master_list)
    return

# create_40_cell_individual_cell_lines()

def create_woc_40_cell():
    save_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/'

    file_paths = ['/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/ccle_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/gcsi_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/gdsc_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/nci60_starfusion_output.csv']

    cell_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/labs/ccle'
    cell_list = [f for f in listdir(cell_dir)]
    csv_list = [['n_datasets_agree', 'cell_line', '#_fusions_predicted', '#_fusions_rejected', 'dataset_list', 'fusion_pred_list', 'fusion_reject_list']]

    for n in range(1, 5):
        for cell in cell_list:
            lab_list = []
            pred_list = []
            pred_list_single = []
            for path in file_paths:
                lab = path.split('/')[-1].split('_')[0]
                with open(path, 'r') as f:
                    reader = csv.reader(f)
                    your_list = list(reader)
                    for line in your_list:
                        if line[1] == cell:
                            pred = ast.literal_eval(line[3])
                            pred_list.append(pred)
                            pred_list_single += pred
                            lab_list.append(lab)

            # get woc
            # get len(pred) before and after to get number of rejections
            # save list before and compare to after (subtract) to get rejected fusions
            pred_list_before_woc = list(dict.fromkeys(pred_list_single))
            pred_list_after_woc = wisdom_of_crowds(pred_list, n)
            num_pred = len(pred_list_after_woc)

            rejected_fusions = list(set(pred_list_before_woc).difference(pred_list_after_woc))
            num_reject = len(rejected_fusions)

            new_list = [n, cell, num_pred, num_reject, lab_list, pred_list_after_woc, rejected_fusions]
            csv_list.append(new_list)

    with open(save_dir + 'woc_starfusion.csv', "w", newline="") as f2:
        writer = csv.writer(f2)
        writer.writerows(csv_list)
    return

# create_woc_40_cell()

def create_confidence_master_list(n_reps_agree=2):

    file_path_40 = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/woc_starfusion.csv'
    file_path_1000 = '/home/m3kowal/Research/genefusions/gene_fusion_master/1000_cell_results/CCLE_Fusions_20181130.txt'
    save_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/1000_cell_results/'

    with open(file_path_1000, 'r') as f:
        reader = csv.reader(f)
        your_list = list(reader)
    fusion_dict_1000 = {}
    for line in your_list[1:]:
        sample_name = line[0].split('\t')[0]
        if sample_name not in fusion_dict_1000.keys():
            fusion_dict_1000[sample_name] = []
            fusion_1000 = line[0].split('\t')[1]
            fusion_dict_1000[sample_name].append(fusion_1000)
        else:
            fusion_1000 = line[0].split('\t')[1]
            fusion_dict_1000[sample_name].append(fusion_1000)

    fusion_list_1000 = []
    for key in fusion_dict_1000:
        no_repeat_list = list(dict.fromkeys(fusion_dict_1000[key]))
        fusion_list_1000 += no_repeat_list

    with open(file_path_40, 'r') as f:
        reader = csv.reader(f)
        your_list = list(reader)
    fusion_list_40 = []

    for line in your_list[1:]:
        if int(line[0]) == n_reps_agree:
            fusion_40 = ast.literal_eval(line[5])
            fusion_list_40 += fusion_40

    fusion_ranking_dict = {}
    for fusion in fusion_list_40:
        fusion_ranking_dict[fusion] = 0

    for fusion in fusion_ranking_dict:
        fusion_count = fusion_list_1000.count(fusion)
        fusion_ranking_dict[fusion] = fusion_count

    ordered_dict = {k: [v, []] for k, v in sorted(fusion_ranking_dict.items(), key=lambda item: item[1], reverse=True)}

    for fusion_40 in ordered_dict:
        for fusion_100 in fusion_dict_1000:
            if fusion_40 in fusion_dict_1000[fusion_100]:
                ordered_dict[fusion_40][1].append(fusion_100)

    csv_fusion_master_list = [['gene_fusion', 'count' , 'cell_lines_found_in']]

    for fusion in ordered_dict:
        new_list = [fusion, ordered_dict[fusion][0], ordered_dict[fusion][1]]
        csv_fusion_master_list.append(new_list)

    with open(save_dir + str(n_reps_agree) + '_reps_agree_vs_1000.csv', "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(csv_fusion_master_list)
    return

# create_confidence_master_list()

def count_occurence_within_40():
    save_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/'

    file_paths = ['/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/ccle_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/gcsi_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/gdsc_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/nci60_starfusion_output.csv']

    cell_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/labs/ccle'
    cell_list = [f for f in listdir(cell_dir)]

    master_fusion_list = []
    fusion_cell_dict = {}
    for cell in cell_list:
        fusion_cell_dict[cell] = []
    lab_dict = {}
    for path in file_paths:
        lab = path.split('/')[7].split('_')[0]
        with open(path, 'r') as f:
            reader = csv.reader(f)
            your_list = list(reader)
            lab_dict[lab] = {}
            for line in your_list[1:]:
                new_list = ast.literal_eval(line[3])
                no_repeat_list = list(dict.fromkeys(new_list))
                lab_dict[lab][line[1]] = no_repeat_list
                master_fusion_list += no_repeat_list

    no_repeat_fusion_master_list = list(dict.fromkeys(master_fusion_list))
    # dict_list = [['fusion', 'occurences', 'cell_lines', 'labs']]
    dict_list = {}
    for fusion in no_repeat_fusion_master_list:
        dict_list[fusion] = [0, [], []]
        for lab in lab_dict:
            for cell in lab_dict[lab]:
                fusion_count = lab_dict[lab][cell].count(fusion)
                if fusion_count > 0:
                    dict_list[fusion][0] += fusion_count
                    dict_list[fusion][1].append(str(lab))
                    dict_list[fusion][2].append(str(cell))


    ordered_dict = {k: v for k, v in sorted(dict_list.items(), key=lambda item: item[1], reverse=True)}
    csv_list = [['fusion', '#_occurences', 'labs_found_in', 'cells_found_in']]
    for fusion in ordered_dict:
        csv_list.append([fusion,
                         ordered_dict[fusion][0],
                         ordered_dict[fusion][1],
                         ordered_dict[fusion][2]])

    with open(save_dir + '40_fusions_count_within.csv', "w", newline="") as f2:
        writer = csv.writer(f2)
        writer.writerows(csv_list)
    return

# count_occurence_within_40()

def count_occurence_within_40_n_agree(n_reps_agree=3):
    save_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/'

    file_paths = ['/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/ccle_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/gcsi_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/gdsc_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/nci60_starfusion_output.csv']

    file_path_40 = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/woc_starfusion.csv'
    cell_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/labs/ccle'
    cell_list = [f for f in listdir(cell_dir)]

    master_fusion_list = []
    fusion_cell_dict = {}
    for cell in cell_list:
        fusion_cell_dict[cell] = []
    lab_dict = {}
    for path in file_paths:
        lab = path.split('/')[7].split('_')[0]
        with open(path, 'r') as f:
            reader = csv.reader(f)
            your_list = list(reader)
            lab_dict[lab] = {}
            for line in your_list[1:]:
                new_list = ast.literal_eval(line[3])
                no_repeat_list = list(dict.fromkeys(new_list))
                lab_dict[lab][line[1]] = no_repeat_list
                master_fusion_list += no_repeat_list

    with open(file_path_40, 'r') as f:
        reader = csv.reader(f)
        your_list = list(reader)
    fusion_list_40 = []

    for line in your_list[1:]:
        if int(line[0]) == n_reps_agree:
            fusion_40 = ast.literal_eval(line[5])
            fusion_list_40 += fusion_40

    no_repeat_fusion_list_40 = list(dict.fromkeys(fusion_list_40))
    # dict_list = [['fusion', 'occurences', 'cell_lines', 'labs']]
    dict_list = {}
    for fusion in no_repeat_fusion_list_40:
        dict_list[fusion] = [0, [], []]
        for lab in lab_dict:
            for cell in lab_dict[lab]:
                fusion_count = lab_dict[lab][cell].count(fusion)
                if fusion_count > 0:
                    dict_list[fusion][0] += fusion_count
                    dict_list[fusion][1].append(str(lab))
                    dict_list[fusion][2].append(str(cell))


    ordered_dict = {k: v for k, v in sorted(dict_list.items(), key=lambda item: item[1], reverse=True)}
    csv_list = [['fusion', '#_occurences', 'labs_found_in', 'cells_found_in']]
    for fusion in ordered_dict:
        csv_list.append([fusion,
                         ordered_dict[fusion][0],
                         ordered_dict[fusion][1],
                         ordered_dict[fusion][2]])

    with open(save_dir + str(n_reps_agree) + '_agree_40_fusions_count_within.csv', "w", newline="") as f2:
        writer = csv.writer(f2)
        writer.writerows(csv_list)
    return

# count_occurence_within_40_n_agree(2)


def occurence_and_evidence_40(n_reps_agree=3, order='JRC'):
    file_path_40 = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/woc_starfusion.csv'
    save_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/evidence/'

    root_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/labs/'

    # iterate through list of 40 (master list)
    with open(file_path_40, 'r') as f:
        reader = csv.reader(f)
        your_list = list(reader)
    fusion_list_40 = []

    for line in your_list[1:]:
        if int(line[0]) == n_reps_agree:
            fusion_40 = ast.literal_eval(line[5])
            fusion_list_40 += fusion_40

    lab_directories = [join(root_dir, f) for f in listdir(root_dir)]


    # create dict = {gene= , evidence1= , evidence2 = , tot_evidence = , occurences = }
    fusion_evidence_dict = {}
    for fusion in fusion_list_40:
        # [ev1, ev2, tot_ev, occurences]
        fusion_evidence_dict[fusion] = [0, 0, 0, [], []]
        for lab_dir in lab_directories:
            sub_dirs = glob(lab_dir + "/*/")
            for subsub_dirs in sub_dirs:
                with open(subsub_dirs + 'star-fusion.fusion_predictions.abridged.tsv') as tsvfile:
                    reader = csv.DictReader(tsvfile, delimiter='\t')
                    for row in reader:
                        if fusion == row['#FusionName']:
                            # if no repeat, add list of current fusions after
                            # 'for row in reader' to keep track of fusions in current doc
                            fusion_evidence_dict[fusion][0] += int(row['JunctionReadCount'])
                            fusion_evidence_dict[fusion][1] += int(row['SpanningFragCount'])
                            fusion_evidence_dict[fusion][2] += int(row['JunctionReadCount'])+int(row['SpanningFragCount'])
    if order == 'JRC':
        ordering_val = 0
    if order == 'SFC':
        ordering_val = 1
    if order == 'TOT':
        ordering_val = 2
    ordered_fusion_evidence_dict = {k: v for k, v in sorted(fusion_evidence_dict.items(), key=lambda item: item[1][ordering_val], reverse=True)}

    file_path_2 = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/40_fusions_count_within.csv'
    with open(file_path_2, 'r') as f:
        reader = csv.reader(f)
        your_list = list(reader)

    for fusion in ordered_fusion_evidence_dict:
        for x in your_list:
            if fusion == x[0]:
                ordered_fusion_evidence_dict[fusion][3] += ast.literal_eval(x[2])
                ordered_fusion_evidence_dict[fusion][4] += ast.literal_eval(x[3])


    csv_list = [['gene', 'JunctionReadCount', 'SpanningFragCount', 'JRC+SFC', 'labs_found_in', 'cells_found_in']]
    for fusion in ordered_fusion_evidence_dict:
        new_list = [fusion, ordered_fusion_evidence_dict[fusion][0],
                    ordered_fusion_evidence_dict[fusion][1],ordered_fusion_evidence_dict[fusion][2],
                    ordered_fusion_evidence_dict[fusion][3],ordered_fusion_evidence_dict[fusion][4]
                    ]
        csv_list.append(new_list)

    with open(save_dir + str(n_reps_agree) + '_agree_40_fusions_JRC_SFC_evidence_order=' + order + '.csv', "w", newline="") as f2:
        writer = csv.writer(f2)
        writer.writerows(csv_list)
    return

# occurence_and_evidence_40(n_reps_agree=3, order='JRC')
# occurence_and_evidence_40(n_reps_agree=3, order='SFC')
# occurence_and_evidence_40(n_reps_agree=3, order='TOT')
#
# occurence_and_evidence_40(n_reps_agree=2, order='JRC')
# occurence_and_evidence_40(n_reps_agree=2, order='SFC')
# occurence_and_evidence_40(n_reps_agree=2, order='TOT')


def count_single_gene_occurence_within_40_n_agree(n_reps_agree=3, partner='3', partner_search='3'):
    save_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/partner_occurence/'

    file_paths = ['/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/ccle_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/gcsi_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/gdsc_starfusion_output.csv',
                  '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/nci60_starfusion_output.csv']

    file_path_40 = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/woc_starfusion.csv'
    cell_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/40_cell_results/labs/ccle'
    cell_list = [f for f in listdir(cell_dir)]

    if partner == '3':
        partner_idx = 0
        other_partner = '5'
    elif partner == '5':
        partner_idx = 1
        other_partner = '3'
    else: exit()

    if partner_search == '3':
        partner_search_idx = 0
    elif partner_search == '5':
        partner_search_idx = 1
    else: exit()

    fusion_cell_dict = {}
    for cell in cell_list:
        fusion_cell_dict[cell] = []
    lab_dict = {}
    for path in file_paths:
        lab = path.split('/')[7].split('_')[0]
        with open(path, 'r') as f:
            reader = csv.reader(f)
            your_list = list(reader)
            lab_dict[lab] = {}
            for line in your_list[1:]:
                new_list = ast.literal_eval(line[3])
                no_repeat_list = list(dict.fromkeys(new_list))
                single_fusion_list = []
                other_fusion_list = []
                for fus in no_repeat_list:
                    single_fusion_list.append(fus.split('--')[0])
                    other_fusion_list.append(fus.split('--')[1])
                lab_dict[lab][line[1]] = {}
                lab_dict[lab][line[1]]['3'] = single_fusion_list
                lab_dict[lab][line[1]]['5'] = other_fusion_list

    with open(file_path_40, 'r') as f:
        reader = csv.reader(f)
        your_list = list(reader)
    fusion_list_40 = []

    for line in your_list[1:]:
        if int(line[0]) == n_reps_agree:
            fusion_40 = ast.literal_eval(line[5])
            single_fusion_40_list = []
            for fus in fusion_40:
                single_fusion_40_list.append(fus.split('--')[partner_idx])
            fusion_list_40 += single_fusion_40_list

    no_repeat_fusion_list_40 = list(dict.fromkeys(fusion_list_40))
    # dict_list = [['fusion', 'occurences', 'cell_lines', 'labs', ]]
    dict_list = {}
    for fusion in no_repeat_fusion_list_40:
        dict_list[fusion] = [0, [], [], []]
        for lab in lab_dict:
            for cell in lab_dict[lab]:
                fusion_count = lab_dict[lab][cell][partner_search].count(fusion)
                if fusion_count > 0:
                    dict_list[fusion][0] += fusion_count
                    dict_list[fusion][1].append(str(lab))
                    dict_list[fusion][2].append(str(cell))
                    i = 0
                    for fus in lab_dict[lab][cell][partner_search]:
                        if fusion == fus:
                            dict_list[fusion][3].append(lab_dict[lab][cell][other_partner][i])
                        i += 1
                    # idx = lab_dict[lab][cell][partner_search].index(fusion)
                    # dict_list[fusion][3].append(lab_dict[lab][cell][other_partner][idx])


    ordered_dict = {k: v for k, v in sorted(dict_list.items(), key=lambda item: item[1], reverse=True)}
    csv_list = [['partner', '#_occurences', 'labs_found_in', 'cells_found_in', 'partners']]

    for fusion in ordered_dict:
        csv_list.append([fusion,
                         ordered_dict[fusion][0],
                         ordered_dict[fusion][1],
                         ordered_dict[fusion][2],
                         ordered_dict[fusion][3],
                         ])

    with open(save_dir + str(n_reps_agree) + '_agree_40_SINGLE_partner=' + str(partner) +
              '_partner_search=' + str(partner_search) + '_fusions_count_within.csv', "w", newline="") as f2:
        writer = csv.writer(f2)
        writer.writerows(csv_list)
    return

count_single_gene_occurence_within_40_n_agree(3, '3', '3')
count_single_gene_occurence_within_40_n_agree(3, '5', '5')
count_single_gene_occurence_within_40_n_agree(3, '3', '5')
count_single_gene_occurence_within_40_n_agree(3, '5', '3')

count_single_gene_occurence_within_40_n_agree(2, '3', '3')
count_single_gene_occurence_within_40_n_agree(2, '5', '5')
count_single_gene_occurence_within_40_n_agree(2, '3', '5')
count_single_gene_occurence_within_40_n_agree(2, '5', '3')