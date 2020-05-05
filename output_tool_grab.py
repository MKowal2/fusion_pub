import sys
import os
import csv
from glob import glob
from evaluations import *


def get_fusioncatcher_pred(path_to_txt):
    pred_tot = []
    pred_first = []
    pred_second = []

    with open(path_to_txt) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            if len(row) > 0:
                if '*' in row[0]:
                    gene_fusion = row[0].split(' ')[3]
                    first_gene = gene_fusion.split('--')[0]
                    second_gene = gene_fusion.split('--')[1]
                    pred_tot.append(gene_fusion)
                    pred_first.append(first_gene)
                    pred_second.append(second_gene)

    pred_lists = [pred_tot, pred_first, pred_second]
    return pred_lists

def get_fusioncatcher_dict(path_to_outputs, output_filename):
    fusion_dict = {}
    directories = glob(path_to_outputs)
    for lab_path in directories:
        subdirectories = glob(lab_path+"/*/")
        labname = lab_path.split('/')[-2]
        fusion_dict[labname] = {}
        for pred_path in subdirectories:
            # get tool predictions in list
            pred_lists = get_fusioncatcher_pred(pred_path + output_filename)
            cell_line = pred_path.split('/')[-2]
            # add predictions to dict
            fusion_dict[labname][cell_line] = pred_lists
    return fusion_dict

# StarFusion - outputs files = star-fusion.fusion_predictions.abridged.tsv / fstar-fusion.fusion_predictions.tsv
def get_tophat_pred(path_to_txt, cell_line):
    pred_tot = []
    pred_first = []
    pred_second = []

    with open(path_to_txt) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            if cell_line == 'BT474-S':
                if 'BT474-S' in row[0] or 'BT747-S' in row[0]:
                    first_gene = row[1]
                    second_gene = row[4]
                    gene_tot = first_gene + '--' + second_gene
                    pred_tot.append(str(gene_tot))
                    pred_first.append(first_gene)
                    pred_second.append(second_gene)
            if cell_line == 'BT474':
                if 'BT474' in row[0] or 'BT747' in row[0]:
                    first_gene = row[1]
                    second_gene = row[4]
                    gene_tot = first_gene + '--' + second_gene
                    pred_tot.append(str(gene_tot))
                    pred_first.append(first_gene)
                    pred_second.append(second_gene)
            else:
                if cell_line in row[0]:
                    first_gene = row[1]
                    second_gene = row[4]
                    gene_tot = first_gene + '--' + second_gene
                    pred_tot.append(str(gene_tot))
                    pred_first.append(first_gene)
                    pred_second.append(second_gene)

    pred_lists = [pred_tot, pred_first, pred_second]
    return pred_lists

def get_tophat_dict(path_to_outputs, output_filename):
    fusion_dict = {}
    directories = glob(path_to_outputs)
    cell_line = ['BT474', 'BT474-S', 'KPL4', 'MCF7', 'SKBR3', 'SKBR3-0', 'SKBR3-1']
    for pred_path in directories:
        labname = pred_path.split('/')[-2]
        fusion_dict[labname] = {}
        # get tool predictions in list
        for cell in cell_line:
            pred_lists = get_tophat_pred(pred_path + output_filename, cell)
            # add predictions to dict
            fusion_dict[labname][cell] = pred_lists
    return fusion_dict

# Arriba - outputs files = fusions.tsv / fusions.discarded.tsv
def get_arriba_pred(path_to_tsv):
    pred_tot = []
    pred_first = []
    pred_second = []

    with open(path_to_tsv) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            first_gene = row['#gene1']
            second_gene = row['gene2']
            if '(' in first_gene:
                if '(' in second_gene:
                    first_split_gene = first_gene.split(',')
                    second_split_gene = second_gene.split(',')
                    first_gene_a = first_split_gene[0].split('(')[0]
                    first_gene_b = first_split_gene[1].split('(')[0]
                    second_gene_a = second_split_gene[0].split('(')[0]
                    second_gene_b = second_split_gene[1].split('(')[0]
                    pred_tot.append(first_gene_a + '--' + second_gene_a)
                    pred_tot.append(first_gene_a + '--' + second_gene_b)
                    pred_tot.append(first_gene_b + '--' + second_gene_a)
                    pred_tot.append(first_gene_b + '--' + second_gene_b)
                    pred_first.append(first_gene_a)
                    pred_first.append(first_gene_a)
                    pred_first.append(first_gene_b)
                    pred_first.append(first_gene_b)
                    pred_second.append(second_gene_a)
                    pred_second.append(second_gene_b)
                    pred_second.append(second_gene_a)
                    pred_second.append(second_gene_b)
                else:
                    first_split_gene = first_gene.split(',')
                    first_gene_a = first_split_gene[0].split('(')[0]
                    first_gene_b = first_split_gene[1].split('(')[0]
                    pred_tot.append(first_gene_a + '--' + second_gene)
                    pred_tot.append(first_gene_b + '--' + second_gene)
                    pred_first.append(first_gene_a)
                    pred_first.append(first_gene_b)
                    pred_second.append(second_gene)
                    pred_second.append(second_gene)
            elif '(' in second_gene:
                second_split_gene = second_gene.split(',')
                second_gene_a = second_split_gene[0].split('(')[0]
                second_gene_b = second_split_gene[1].split('(')[0]
                pred_tot.append(first_gene + '--' + second_gene_a)
                pred_tot.append(first_gene + '--' + second_gene_b)
                pred_first.append(first_gene)
                pred_first.append(first_gene)
                pred_second.append(second_gene_a)
                pred_second.append(second_gene_b)
            else:
                pred_tot.append(first_gene + '--' + second_gene)
                pred_first.append(first_gene)
                pred_second.append(second_gene)

    pred_lists = [pred_tot, pred_first, pred_second]
    return pred_lists

def get_arriba_dict(path_to_outputs, output_filename):
    fusion_dict = {}
    directories = glob(path_to_outputs)
    for lab_path in directories:
        subdirectories = glob(lab_path+"/*/")
        labname = lab_path.split('/')[-2]
        fusion_dict[labname] = {}
        for pred_path in subdirectories:
            # get tool predictions in list
            pred_lists = get_arriba_pred(pred_path + output_filename)
            cell_line = pred_path.split('/')[-2]
            # add predictions to dict
            fusion_dict[labname][cell_line] = pred_lists
    return fusion_dict

# StarFusion - outputs files = star-fusion.fusion_predictions.abridged.tsv / fstar-fusion.fusion_predictions.tsv

def get_starfusion_pred(path_to_tsv):
    pred_tot = []
    pred_first = []
    pred_second = []

    with open(path_to_tsv) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            fusion_pred = row['#FusionName']
            pred_tot.append(str(fusion_pred))
            pred_first.append(fusion_pred.split('--')[0])
            pred_second.append(fusion_pred.split('--')[1])

    pred_lists = [pred_tot, pred_first, pred_second]
    return pred_lists

def get_starfusion_dict(path_to_outputs, output_filename):
    fusion_dict = {}
    directories = glob(path_to_outputs)
    for lab_path in directories:
        subdirectories = glob(lab_path+"/*/")
        labname = lab_path.split('/')[-2]
        fusion_dict[labname] = {}
        for pred_path in subdirectories:
            # get tool predictions in list
            pred_lists = get_starfusion_pred(pred_path + output_filename)
            cell_line = pred_path.split('/')[-2]
            # add predictions to dict
            fusion_dict[labname][cell_line] = pred_lists
    return fusion_dict
