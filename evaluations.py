import sys
import os
import csv
from glob import glob


def return_gt_fusions(path):
    '''
    input: xlrd.open_workbook('path_to_workbook.xlsx')
    returns: fusion list ['gene1--gene2', ...], first gene list ['gene1', ...], second gene list ['gene2', ...],
    '''
    # assert cell_line in ['BT474', 'KPL4', 'MCF7', 'SKBR3'], 'Cell line not one of the four experimentally valitated!'

    # open workbook and grab sheet
    # sheet = workbook.sheet_by_name('Table S3 - ExpValidFusions')
    fusion_gt_dict = {}
    # bt474
    with open(path + 'bt474_gt.csv') as file:
        bt474_list = [line.rstrip('\n').split(',')[0] for line in file]
    fusion_gt_dict['BT474'] = bt474_list

    # kpl4
    with open(path + 'kpl4_gt.csv') as file:
        kpl4_list = [line.rstrip('\n').split(',')[0] for line in file]
    fusion_gt_dict['KPL4'] = kpl4_list

    # mcf7
    with open(path + 'mcf7_gt.csv') as file:
        mcf7_list = [line.rstrip('\n').split(',')[0] for line in file]
    fusion_gt_dict['MCF7'] = mcf7_list

    # skbr3
    with open(path + 'skbr3_gt.csv') as file:
        skbr3_list = [line.rstrip('\n').split(',')[0] for line in file]
    fusion_gt_dict['SKBR3'] = skbr3_list

    return fusion_gt_dict

def return_gt_fusions_copy(workbook):
    '''
    input: xlrd.open_workbook('path_to_workbook.xlsx')
    returns: fusion list ['gene1--gene2', ...], first gene list ['gene1', ...], second gene list ['gene2', ...],
    '''
    assert cell_line in ['BT474', 'KPL4', 'MCF7', 'SKBR3'], 'Cell line not one of the four experimentally valitated!'

    # open workbook and grab sheet
    # sheet = workbook.sheet_by_name('Table S3 - ExpValidFusions')

    # start of fusion list
    row_start = 107
    row_end = 160
    first_gene_list = []
    second_gene_list = []
    # seperated by --
    fusion_list = []
    for i in range(row_start, row_end):
        fusion = sheet.cell(i, 1).value
        if cell_line in fusion:
            fusion = fusion[len(cell_line)+1:]
            fusion_list.append(fusion)
            first_gene_list.append(fusion.split('--')[0])
            second_gene_list.append(fusion.split('--')[1])
    return fusion_list, first_gene_list, second_gene_list


# EVALUATION tool predictions against gt
def evaluate_tool(pred, cell_line, gt, print_fusion_list=False):
    '''
    inputs - pred = list of predicted fusions
    cell_line = string of cell line
    gt - list of gene fusions ['gene1--gene2', ...]
    returns - true positives, accuracy and list of fusions predicted correctly
    '''
    assert cell_line in ['BT474', 'BT474-S', 'KPL4', 'MCF7', 'SKBR3', 'SKBR3-0', 'SKBR3-1'], 'Cell line not one of the four experimentall valitated!'

    # gt = [bt474_gt_fusions, kpl4_gt_fusions, mcf7_gt_fusions, skbr3_gt_fusions]

    if cell_line == 'BT474':
        gt = gt['BT474']
    elif cell_line == 'BT474-S':
        gt = gt['BT474']
    elif cell_line == 'KPL4':
        gt = gt['KPL4']
    elif cell_line == 'MCF7':
        gt = gt['MCF7']
    elif cell_line == 'SKBR3':
        gt = gt['SKBR3']
    elif cell_line == 'SKBR3-0':
        gt = gt['SKBR3']
    elif cell_line == 'SKBR3-1':
        gt = gt['SKBR3']

    pred_combined = pred
    gt_combined = gt

    true_positives = []
    false_positives = []
    for tool_pred in pred_combined:
        if tool_pred in gt_combined:
            true_positives.append(tool_pred)
        else:
            false_positives.append(tool_pred)


    true_positives = list(dict.fromkeys(true_positives))
    print(true_positives)
    false_positives = list(dict.fromkeys(false_positives))
    accuracy = len(true_positives)/len(gt_combined)
    if (len(true_positives)+len(false_positives)) > 0:
        precision = len(true_positives)/ (len(true_positives)+len(false_positives))
    else:
        precision = 0
    print('----------------------------------------')
    print('The accuracy is {}/{} = {:.4f}%'.format(len(true_positives), len(gt_combined), accuracy*100))
    print('The precision is {}/({}+{}) = {:.4f}%'.format(len(true_positives),
                                                    len(true_positives),len(false_positives), precision*100))
    # if print_fusion_list:
    #     print('The tool guessed these fusions correctly:')
    #     print(true_positives)
    print('----------------------------------------')

    return [cell_line, accuracy, precision, len(true_positives), len(gt_combined), len(false_positives), true_positives, false_positives]


def wisdom_of_crowds(fusion_lists, n_tools=1):
    '''
    inputs:
    fusion_lists = list of lists (each element is gene_fusion) //
    - n_tools = number of tools to agree
    returns:
    - woc_output = list of all fusions with  >= n_tools agreements
    '''

    # create total list -> remove dups -> iterate through list and check for n agreements -> append to list of fusions
    fusion_master_list = []
    for fusion_list in fusion_lists:
        fusion_master_list += fusion_list

    fusion_master_list = list(dict.fromkeys(fusion_master_list))

    woc_output = []
    for fusion in fusion_master_list:
        n_agree = 0
        for fusion_list in fusion_lists:
            if fusion in fusion_list:
                n_agree += 1
        if n_agree > n_tools-1:
            woc_output.append(fusion)
    return woc_output
