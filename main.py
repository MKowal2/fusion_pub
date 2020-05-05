import csv
from itertools import combinations
from output_tool_grab import *
from evaluations import *

######################################################################################################################
######################################################################################################################

def evaluate_single_tool_single_dataset_write_csv(pred, tool):
    '''
    Inputs:
    pred - tool prediction dictionary
    tool - name of tool being evaluated
    :return:
    outputs csv file of results to: './outputs/1-tool-1-dataset/' + tool + "_individual_datasets.csv"
    '''

    # take output list of tool for that cell line / dataset
    # evaluate against gt list
    # grab ground truth fusions for all cell lines
    gt = return_gt_fusions('./validated_fusions/')

    csv_list = [['dataset', 'cell', 'accuracy', 'precision', '#TP', '#GT', '#FP', 'TP_list', 'FP_list']]

    #ccle
    lab = 'ccle'
    eval_list = evaluate_tool(pred[lab]['BT474'][0], 'BT474', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    eval_list = evaluate_tool(pred[lab]['MCF7'][0], 'MCF7', gt,  print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    eval_list = evaluate_tool(pred[lab]['SKBR3'][0], 'SKBR3', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    #edgren
    lab = 'edgren'
    eval_list = evaluate_tool(pred[lab]['BT474'][0], 'BT474', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    eval_list = evaluate_tool(pred[lab]['BT474-S'][0], 'BT474-S', gt, print_fusion_list=True)
    eval_list.insert(0, lab)
    csv_list.append(eval_list)

    eval_list = evaluate_tool(pred[lab]['MCF7'][0], 'MCF7', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    eval_list = evaluate_tool(pred[lab]['SKBR3-0'][0], 'SKBR3-0', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    eval_list = evaluate_tool(pred[lab]['SKBR3-1'][0], 'SKBR3-1', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    # eval_list = evaluate_tool(pred[lab]['KPL4'][0], 'KPL4', gt, print_fusion_list=True)
    # eval_list.insert(0,lab)
    # csv_list.append(eval_list)

    if tool == 'tophat':
        #edgren
        lab = 'edgren_lenient'
        eval_list = evaluate_tool(pred[lab]['BT474'][0], 'BT474', gt, print_fusion_list=True)
        eval_list.insert(0,lab)
        csv_list.append(eval_list)

        # eval_list = evaluate_tool(pred[lab]['BT474-S'][0], 'BT474-S', gt, print_fusion_list=True)
        # eval_list.insert(0, lab)
        # csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['MCF7'][0], 'MCF7', gt, print_fusion_list=True)
        eval_list.insert(0,lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['SKBR3-0'][0], 'SKBR3-0', gt, print_fusion_list=True)
        eval_list.insert(0,lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['SKBR3-1'][0], 'SKBR3-1', gt, print_fusion_list=True)
        eval_list.insert(0,lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['KPL4'][0], 'KPL4', gt, print_fusion_list=True)
        eval_list.insert(0,lab)
        csv_list.append(eval_list)

    #gcsi
    lab = 'gcsi'
    eval_list = evaluate_tool(pred[lab]['BT474'][0], 'BT474', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    eval_list = evaluate_tool(pred[lab]['MCF7'][0], 'MCF7', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    eval_list = evaluate_tool(pred[lab]['SKBR3'][0], 'SKBR3', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    eval_list = evaluate_tool(pred[lab]['KPL4'][0], 'KPL4', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    if tool == 'tophat':
        # gcsi lenient
        lab = 'gcsi_lenient'
        eval_list = evaluate_tool(pred[lab]['BT474'][0], 'BT474', gt, print_fusion_list=True)
        eval_list.insert(0, lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['MCF7'][0], 'MCF7', gt, print_fusion_list=True)
        eval_list.insert(0, lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['SKBR3'][0], 'SKBR3', gt, print_fusion_list=True)
        eval_list.insert(0, lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['KPL4'][0], 'KPL4', gt, print_fusion_list=True)
        eval_list.insert(0, lab)
        csv_list.append(eval_list)

    # if tool == 'tophat':
        # ccle lenient
        # lab = 'ccle_lenient'
        # eval_list = evaluate_tool(pred[lab]['BT474'][0], 'BT474', gt, print_fusion_list=True)
        # eval_list.insert(0, lab)
        # csv_list.append(eval_list)
        #
        # eval_list = evaluate_tool(pred[lab]['MCF7'][0], 'MCF7', gt, print_fusion_list=True)
        # eval_list.insert(0, lab)
        # csv_list.append(eval_list)
        #
        # eval_list = evaluate_tool(pred[lab]['SKBR3'][0], 'SKBR3', gt, print_fusion_list=True)
        # eval_list.insert(0, lab)
        # csv_list.append(eval_list)

    #gray
    if tool == 'tophat':
        lab = 'gray'
        eval_list = evaluate_tool(pred[lab]['BT474'][0], 'BT474', gt, print_fusion_list=True)
        eval_list.insert(0,lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['MCF7'][0], 'MCF7', gt, print_fusion_list=True)
        eval_list.insert(0,lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['SKBR3'][0], 'SKBR3', gt, print_fusion_list=True)
        eval_list.insert(0,lab)
        csv_list.append(eval_list)
        # grey lenient
        lab = 'gray_lenient'
        eval_list = evaluate_tool(pred[lab]['BT474'][0], 'BT474', gt, print_fusion_list=True)
        eval_list.insert(0, lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['MCF7'][0], 'MCF7', gt, print_fusion_list=True)
        eval_list.insert(0, lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['SKBR3'][0], 'SKBR3', gt, print_fusion_list=True)
        eval_list.insert(0, lab)
        csv_list.append(eval_list)

    else:
        lab = 'gray'
        eval_list = evaluate_tool(pred[lab]['BT474'][0], 'BT474', gt, print_fusion_list=True)
        eval_list.insert(0, lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['MCF7'][0], 'MCF7', gt, print_fusion_list=True)
        eval_list.insert(0, lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['SKBR3'][0], 'SKBR3', gt, print_fusion_list=True)
        eval_list.insert(0, lab)
        csv_list.append(eval_list)

    #uhn
    lab = 'uhn'
    eval_list = evaluate_tool(pred[lab]['BT474'][0], 'BT474', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    eval_list = evaluate_tool(pred[lab]['MCF7'][0], 'MCF7', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    eval_list = evaluate_tool(pred[lab]['SKBR3'][0], 'SKBR3', gt, print_fusion_list=True)
    eval_list.insert(0,lab)
    csv_list.append(eval_list)

    if tool == 'tophat':
        #uhn
        lab = 'uhn_lenient'
        eval_list = evaluate_tool(pred[lab]['BT474'][0], 'BT474', gt, print_fusion_list=True)
        eval_list.insert(0,lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['MCF7'][0], 'MCF7', gt, print_fusion_list=True)
        eval_list.insert(0,lab)
        csv_list.append(eval_list)

        eval_list = evaluate_tool(pred[lab]['SKBR3'][0], 'SKBR3', gt, print_fusion_list=True)
        eval_list.insert(0,lab)
        csv_list.append(eval_list)

    with open('./outputs/1-tool-1-dataset/' + tool + "_individual_datasets.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(csv_list)
    return

######################################################################################################################
######################################################################################################################

def evaluate_woc_tool_single_dataset_write_csv(pred_dict, tool_list):
    '''
    Inputs:
    pred - tool prediction dictionary
    tool_list - list of tools being evaluated
    :return:
    outputs csv file of results to: './outputs/woc_tools_1-dataset/woc_tools_individual_datasets_V2.csv'
    '''

    # grab ground truth fusions for all cell lines
    gt = return_gt_fusions('./validated_fusions/')

    csv_list = [['ntools','dataset', 'cell', 'accuracy', 'precision', '#TP', '#GT', '#FP', 'TP_list', 'FP_list']]
    labs = {}
    labs['ccle'] = ['BT474', 'MCF7', 'SKBR3']
    # labs['edgren'] = ['BT474', 'BT474-S', 'KPL4', 'MCF7', 'SKBR3-0', 'SKBR3-1']
    labs['edgren'] = ['BT474', 'KPL4', 'MCF7', 'SKBR3']
    labs['gcsi'] = ['BT474', 'KPL4', 'MCF7', 'SKBR3']
    labs['gray'] = ['BT474', 'MCF7', 'SKBR3', ]
    labs['uhn'] = ['BT474', 'MCF7', 'SKBR3']

    for n in range(1, len(tool_list)+1):
        for lab in labs:
            for cell in labs[lab]:
                pred = []
                for tool in pred_dict:
                    if lab == 'edgren' and cell == 'SKBR3':
                        single_pred = pred_dict[tool][lab]['SKBR3-0'][0]
                        pred.append(single_pred)
                    else:
                        single_pred = pred_dict[tool][lab][cell][0]
                        pred.append(single_pred)

                woc_output = wisdom_of_crowds(pred, n)
                eval_list = evaluate_tool(woc_output, cell, gt)
                eval_list.insert(0, lab)
                eval_list.insert(0, n)
                csv_list.append(eval_list)

    csv_list.insert(0, [])
    csv_list.insert(0, tool_list)
    with open('./outputs/woc_tools_1-dataset/woc_tools_individual_datasets_V2.csv', "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(csv_list)
    return

######################################################################################################################
######################################################################################################################

def evaluate_single_tool_woc_dataset_write_csv(pred_dict, tool, lab_list, sk_val_to_use,use_short_bt=False):
    '''
    Inputs:
    pred - tool prediction dictionary
    tool - tool being evaluated
    lab_list - list of labs being evaluated
    :return:
    outputs csv file of results to: './outputs/1-tool_woc_dataset/' + tool + '_skb' + str(sk_val_to_use) + '_1_tools_woc_datasets.csv'
    (for case of not using SKBR3-1 and not using BT474-S)
    '''

    # grab ground truth fusions for all cell lines
    gt = return_gt_fusions('./validated_fusions/')

    csv_list = [['ndatasets', 'cell', 'accuracy', 'precision', '#TP', '#GT', '#FP', 'TP_list', 'FP_list']]

    cells = {}
    cells['BT474'] = ['ccle', 'gcsi', 'gray', 'edgren', 'uhn']
    # cells['BT474-S'] = ['edgren']
    cells['KPL4'] = ['gcsi', 'edgren']
    cells['MCF7'] = ['ccle', 'gcsi', 'gray', 'edgren', 'uhn']
    cells['SKBR3'] = ['ccle', 'gcsi', 'edgren', 'gray', 'uhn']
    # cells['SKBR3-0'] = ['edgren']
    # cells['SKBR3-1'] = ['edgren']

    for n in range(1, len(lab_list)+1):
        for cell in cells:
            pred = []
            labs_that_agree = []
            for lab in cells[cell]:
                if lab == 'edgren' and cell == 'BT474':
                    if use_short_bt:
                        single_pred = pred_dict[tool][lab]['BT474-S'][0]
                        pred.append(single_pred)
                    else:
                        single_pred = pred_dict[tool][lab][cell][0]
                        pred.append(single_pred)
                elif lab == 'edgren' and cell == 'SKBR3':
                        if sk_val_to_use == 0:
                            single_pred = pred_dict[tool][lab]['SKBR3-0'][0]
                            pred.append(single_pred)
                        elif sk_val_to_use == 1:
                            single_pred = pred_dict[tool][lab]['SKBR3-1'][0]
                            pred.append(single_pred)
                else:
                    single_pred = pred_dict[tool][lab][cell][0]
                    pred.append(single_pred)

            woc_output = wisdom_of_crowds(pred, n)
            eval_list = evaluate_tool(woc_output, cell, gt)
            eval_list.insert(0, n)
            csv_list.append(eval_list)

    csv_list.insert(0, [])
    csv_list.insert(0, lab_list)
    if use_short_bt:
        with open('./outputs/1-tool_woc_dataset/' + tool + '_skb' + str(sk_val_to_use) + '_bt-s_1_tools_woc_datasets.csv', "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(csv_list)
    else:
        with open(
                './outputs/1-tool_woc_dataset/' + tool + '_skb' + str(sk_val_to_use) + '_1_tools_woc_datasets.csv',
                "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(csv_list)

    return

######################################################################################################################
######################################################################################################################

def evaluate_woc_all_tool_woc_all_dataset_nagree_write_csv(pred_dict, tool_list, lab_list, sk_val_to_use=1, use_short_bt=False, bytool=False):

    gt = return_gt_fusions('./validated_fusions/')

    # which tools used / wich labs used
    csv_list = [['ndatasets_agree', 'ntools_agree', 'cell', 'accuracy', 'precision', '#TP', '#GT', '#FP', 'TP_list', 'FP_list']]

    cells = {}
    cells['BT474'] = ['ccle', 'gcsi', 'gray', 'edgren', 'uhn']
    # cells['BT474-S'] = ['edgren']
    cells['KPL4'] = ['gcsi', 'edgren']
    cells['MCF7'] = ['ccle', 'gcsi', 'gray', 'edgren', 'uhn']
    cells['SKBR3'] = ['ccle', 'gcsi', 'edgren', 'gray', 'uhn']

    # create master dictionaries for tools, datasets, and tools-pre wisdom of crowds
    tool_master_dataset_dict = {}
    woc_master_dataset_dict = {}
    pre_woc_tool = {}
    for cell in cells:
        tool_master_dataset_dict[cell] = {}
        woc_master_dataset_dict[cell] = {}
        pre_woc_tool[cell] = {}
        for tool_name in tool_list:
            tool_master_dataset_dict[cell][tool_name] = []
            woc_master_dataset_dict[cell][tool_name] = []

    # create master dict per tool (across labs)
    for tool in tool_list:
        for cell in cells:
            # get master list of woc for n_datasets for a single cell
            for lab in cells[cell]:
                # taking care of outlier cases
                if lab == 'edgren' and cell == 'BT474':
                    if use_short_bt:
                        single_pred = pred_dict[tool][lab]['BT474-S'][0]
                        tool_master_dataset_dict[cell][tool].append(single_pred)
                    else:
                        single_pred = pred_dict[tool][lab][cell][0]
                        tool_master_dataset_dict[cell][tool].append(single_pred)
                # taking care of outlier cases
                elif lab == 'edgren' and cell == 'SKBR3':
                    if sk_val_to_use == 0:
                        single_pred = pred_dict[tool][lab]['SKBR3-0'][0]
                        tool_master_dataset_dict[cell][tool].append(single_pred)
                    elif sk_val_to_use == 1:
                        single_pred = pred_dict[tool][lab]['SKBR3-1'][0]
                        tool_master_dataset_dict[cell][tool].append(single_pred)
                else:
                    single_pred = pred_dict[tool][lab][cell][0]
                    tool_master_dataset_dict[cell][tool].append(single_pred)

    # do wisdom of crowds across datasets
    for cell in cells:
        for tool in tool_list:
            for n_data in range(1, len(lab_list)+1):
                woc_master_dataset_dict[cell][tool].append(wisdom_of_crowds(tool_master_dataset_dict[cell][tool], n_data))

    # to wisdom of crowds across tools
    for cell in cells:
        for n in range(len(lab_list)):
            pre_woc_tool[cell][n] = []
            for tool in tool_list:
                pre_woc_tool[cell][n].append(woc_master_dataset_dict[cell][tool][n])

    #
    woc_master_list = {}
    for cell in cells:
        woc_master_list[cell] = [[[],[],[],[],[]] , [[],[],[],[],[]], [[],[],[],[],[]], [[],[],[],[],[]]]

    for cell in cells:
        for n_datasets in range(1, len(pre_woc_tool[cell])+1):
            for n_tools in range(1,len(tool_list)+1):
                woc_master_list[cell][n_tools-1][n_datasets-1] = wisdom_of_crowds(pre_woc_tool[cell][n_datasets-1], n_tools)


    for cell in cells:
        # output file in order of tools 1-n or by labs 1-m
        if bytool:
            for n_datasets in range(len(lab_list)):
                for n_tools in range(len(tool_list)):
                    pred_list = evaluate_tool(woc_master_list[cell][n_tools][n_datasets], cell, gt)
                    pred_list.insert(0, n_tools+1)
                    pred_list.insert(0, n_datasets+1)
                    csv_list.append(pred_list)
        else:
            for n_tools in range(len(tool_list)):
                for n_datasets in range(len(lab_list)):
                    pred_list = evaluate_tool(woc_master_list[cell][n_tools][n_datasets], cell, gt)
                    pred_list.insert(0, n_tools+1)
                    pred_list.insert(0, n_datasets+1)
                    csv_list.append(pred_list)

    # insert labs used and tools used
    csv_list.insert(0, [])
    csv_list.insert(0, lab_list)
    csv_list.insert(0, tool_list)

    # save csv
    if use_short_bt:
        with open('./outputs/woc-tool_woc-dataset/' + 'skb' + str(sk_val_to_use) + '_bt-s_woc_tools_woc_datasets_BYTOOLS_V2.csv', "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(csv_list)
    else:
        with open('./outputs/woc-tool_woc-dataset/' + 'skb' + str(sk_val_to_use) + '_woc_tools_woc_datasets_BYTOOLS_V2.csv', "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(csv_list)
    return

######################################################################################################################
######################################################################################################################

def withold_single_tool_woc_dataset_write_csv(pred_dict, tool, lab_list, sk_val_to_use,use_short_bt=False):
    '''
    Inputs:
     - output of tool
     - name of tool being evaluated
    :return:
    - writes csv file
    '''

    # grab ground truth fusions for all cell lines
    gt = return_gt_fusions('./validated_fusions/')

    csv_list = [['ndatasets', 'cell', 'accuracy', 'precision', '#TP', '#GT', '#FP', 'TP_list', 'FP_list']]
    # remove unused labs from list
    lab_list_removal = ['ccle', 'gcsi', 'gray', 'edgren', 'uhn']
    for lab in lab_list:
        lab_list_removal.remove(lab)

    cells = {}
    cells['BT474'] = ['ccle', 'gcsi', 'gray', 'edgren', 'uhn']
    # cells['BT474-S'] = ['edgren']
    # cells['KPL4'] = ['gcsi', 'edgren']
    cells['MCF7'] = ['ccle', 'gcsi', 'gray', 'edgren', 'uhn']
    cells['SKBR3'] = ['ccle', 'gcsi', 'edgren', 'gray', 'uhn']
    # cells['SKBR3-0'] = ['edgren']
    # cells['SKBR3-1'] = ['edgren']

    # remove unused labs from cell dictionary
    for remove_lab in lab_list_removal:
        for cell in cells:
            if remove_lab in cells[cell]: cells[cell].remove(remove_lab)

    for n in range(1, len(lab_list)+1):
        for cell in cells:
            pred = []
            for lab in cells[cell]:
                # taking care of outlier cases
                if lab == 'edgren' and cell == 'BT474':
                    if use_short_bt:
                        single_pred = pred_dict[tool][lab]['BT474-S'][0]
                        pred.append(single_pred)
                    else:
                        single_pred = pred_dict[tool][lab][cell][0]
                        pred.append(single_pred)
                # taking care of outlier cases
                elif lab == 'edgren' and cell == 'SKBR3':
                        if sk_val_to_use == 0:
                            single_pred = pred_dict[tool][lab]['SKBR3-0'][0]
                            pred.append(single_pred)
                        elif sk_val_to_use == 1:
                            single_pred = pred_dict[tool][lab]['SKBR3-1'][0]
                            pred.append(single_pred)
                # taking care of NORMAL cases
                else:
                    single_pred = pred_dict[tool][lab][cell][0]
                    pred.append(single_pred)
            # add results to overall list
            woc_output = wisdom_of_crowds(pred, n)
            eval_list = evaluate_tool(woc_output, cell, gt)
            eval_list.insert(0, n)
            csv_list.append(eval_list)
    # insert blank line and labs used to overall list
    csv_list.insert(0, [])
    csv_list.insert(0, lab_list)
    # save csv file
    if use_short_bt:
        with open('./outputs/withold_1-tool_woc_dataset/' + str(len(lab_list)) + '_datasets_agree/'
                  + tool + str(lab_list) + '_skb' + str(sk_val_to_use) + '_bt-s_1_tools_woc_datasets.csv', "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(csv_list)
    else:
        with open('./outputs/withold_1-tool_woc_dataset/' + str(len(lab_list)) + '_datasets_agree/'
                  + tool + str(lab_list) + '_skb' + str(sk_val_to_use) + '_1_tools_woc_datasets.csv', "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(csv_list)
    return

######################################################################################################################
######################################################################################################################

def withold_woc_tool_single_dataset_write_csv(pred_dict, lab, tool_list, sk_val_to_use,use_short_bt=False):
    '''
    Inputs:
     - output of tool
     - name of tool being evaluated
    :return:
    - writes csv file
    '''

    # grab ground truth fusions for all cell lines
    gt = return_gt_fusions('./validated_fusions/')

    csv_list = [['ntools', 'lab', 'cell', 'accuracy', 'precision', '#TP', '#GT', '#FP', 'TP_list', 'FP_list']]

    labs = {}
    labs['ccle'] = ['BT474', 'MCF7', 'SKBR3']
    # labs['edgren'] = ['BT474', 'BT474-S', 'KPL4', 'MCF7', 'SKBR3-0', 'SKBR3-1']
    labs['edgren'] = ['BT474', 'KPL4', 'MCF7', 'SKBR3']
    labs['gcsi'] = ['BT474', 'KPL4', 'MCF7', 'SKBR3']
    labs['gray'] = ['BT474', 'MCF7', 'SKBR3', ]
    labs['uhn'] = ['BT474', 'MCF7', 'SKBR3']


    for n in range(1, len(tool_list)+1):
        for cell in labs[lab]:
            pred = []
            tools_that_agree = []
            for tool in tool_list:
                # taking care of outlier cases
                if lab == 'edgren' and cell == 'BT474':
                    if use_short_bt:
                        single_pred = pred_dict[tool][lab]['BT474-S'][0]
                        pred.append(single_pred)
                    else:
                        single_pred = pred_dict[tool][lab][cell][0]
                        pred.append(single_pred)
                # taking care of outlier cases
                elif lab == 'edgren' and cell == 'SKBR3':
                        if sk_val_to_use == 0:
                            single_pred = pred_dict[tool][lab]['SKBR3-0'][0]
                            pred.append(single_pred)
                        elif sk_val_to_use == 1:
                            single_pred = pred_dict[tool][lab]['SKBR3-1'][0]
                            pred.append(single_pred)
                # taking care of NORMAL cases
                else:
                    single_pred = pred_dict[tool][lab][cell][0]
                    pred.append(single_pred)
                # add results to overall list
            woc_output = wisdom_of_crowds(pred, n)
            eval_list = evaluate_tool(woc_output, cell, gt)
            eval_list.insert(0, lab)
            eval_list.insert(0, n)
            csv_list.append(eval_list)
    # insert blank line and labs used to overall list
    csv_list.insert(0, [])
    csv_list.insert(0, tool_list)
    # save csv file
    if use_short_bt:
        with open('./outputs/withold_woc_tool_1-dataset/' + str(len(tool_list)) + '_tools_agree/'
                  + lab + str(tool_list) + '_skb' + str(sk_val_to_use) + '_bt-s_1_tools_woc_datasets.csv', "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(csv_list)
    else:
        with open('./outputs/withold_woc_tool_1-dataset/' + str(len(tool_list)) + '_tools_agree/'
                  + lab + str(tool_list) + '_skb' + str(sk_val_to_use) + '_1_tools_woc_datasets.csv', "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(csv_list)
    return

######################################################################################################################
######################################################################################################################

# obtain predictions
fusioncatcher_path = './tool_outputs_WOC/Fusioncatcher/*/'
fusioncatcher_pred = get_fusioncatcher_dict(fusioncatcher_path, 'summary_candidate_fusions.txt')

starfusion_path = './tool_outputs_WOC/Star_fusion/*/'
starfusion_pred = get_starfusion_dict(starfusion_path, 'star-fusion.fusion_predictions.abridged.tsv')

tophat_path = './tool_outputs_WOC/tophat/*/'
tophat_pred = get_tophat_dict(tophat_path, 'result.txt')

arriba_path = './tool_outputs_WOC/arriba/*/'
arriba_pred = get_arriba_dict(arriba_path, 'fusions.tsv')

# create pre_dict from each tools predictions
pred_dict = {}
pred_dict['fusioncatcher'] = fusioncatcher_pred
pred_dict['starfusion'] = starfusion_pred
pred_dict['tophat'] = tophat_pred
pred_dict['arriba'] = arriba_pred
# define tool list / lab list
tool_list = ['fusioncatcher', 'starfusion', 'tophat', 'arriba']
lab_list = ['ccle', 'gcsi', 'gray', 'edgren', 'uhn']

######################################################################################################################
######################################################################################################################


# running all 4 functions // uncomment to run!
#######################################################################################################################
# 1 tool - 1 dataset
# evaluate_single_tool_single_dataset_write_csv(fusioncatcher_pred, tool='fusioncatcher')
# evaluate_single_tool_single_dataset_write_csv(arriba_pred, tool='arriba')
# evaluate_single_tool_single_dataset_write_csv(tophat_pred, tool='tophat')
# evaluate_single_tool_single_dataset_write_csv(starfusion_pred, tool='starfusion')
#
# # woc tool - 1 dataset
# evaluate_woc_tool_single_dataset_write_csv(pred_dict, tool_list)
#
#
# # 1 tool - woc datasets
# for tool in tool_list:
#     evaluate_single_tool_woc_dataset_write_csv(pred_dict, tool, lab_list, sk_val_to_use=0, use_short_bt=False)
#
# # woc tool - woc dataset
# evaluate_woc_all_tool_woc_all_dataset_nagree_write_csv(pred_dict, tool_list, lab_list, sk_val_to_use = 0,use_short_bt=False)
#
# # withold - 1 tool - woc dataset
# # generate all combinations of labs (n=3)
# permutation_lab_list = list(combinations(lab_list, 4))
# sk_val_to_use = 0
# for tool in tool_list:
#     for lab_list in permutation_lab_list:
#         withold_single_tool_woc_dataset_write_csv(pred_dict, tool, lab_list, sk_val_to_use, use_short_bt=False)


permutation_tool_list = list(combinations(tool_list, 4))
sk_val_to_use = 0
for lab in lab_list:
    for list_of_tools in permutation_tool_list:
        withold_woc_tool_single_dataset_write_csv(pred_dict, lab, list_of_tools, sk_val_to_use, use_short_bt=False)
######################################################################################################################