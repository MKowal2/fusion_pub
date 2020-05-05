import sys
import os
import csv
from os import listdir
from os.path import isfile, join
from glob import glob
import numpy as np
from itertools import combinations
from output_tool_grab import *
from evaluations import *

def create_cell_average_withold():
    path = '/home/m3kowal/Research/genefusions/gene_fusion_master/outputs/withold_1-tool_woc_dataset/2_datasets_agree/'
    save_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/outputs/withold_1-tool_woc_dataset/2_datasets_agree/cell_line_averages/'


    files = [f for f in listdir(path) if isfile(join(path, f))]
    for file in files:
        with open(path + file, 'r') as f:
            reader = csv.reader(f)
            your_list = list(reader)
            csv_list = [['ndatasets', 'accuracy', 'precision', '#Total-TP', '#Total-GT', '#Total-FP']]
            for n in range(1, 3):
                new_csv_list = [0.0, 0.0, 0, 0, 0]
                for line_list in your_list:
                    if len(line_list) > 0:
                        if line_list[0] == str(n):
                            new_csv_list[0] += float(line_list[2])
                            new_csv_list[1] += float(line_list[3])
                            new_csv_list[2] += int(line_list[4])
                            new_csv_list[3] += int(line_list[5])
                            new_csv_list[4] += int(line_list[6])
                new_csv_list[0] = new_csv_list[0]/3
                new_csv_list[1] = new_csv_list[1]/3
                new_csv_list.insert(0, n)
                csv_list.append(new_csv_list)

        with open(save_dir + file, "w", newline="") as f2:
            writer = csv.writer(f2)
            writer.writerows(csv_list)
    return


def create_cell_average_woc_tool_indiv_dataset():
    csv_file = '/home/m3kowal/Research/genefusions/gene_fusion_master/outputs/woc_tools_1-dataset/woc_tools_individual_datasets_V2.csv'
    save_dir = '/home/m3kowal/Research/genefusions/gene_fusion_master/outputs/woc_tools_1-dataset/'
    csv_list = [['ntools, dataset', 'accuracy', 'precision', '#Total-TP', '#Total-GT', '#Total-FP']]
    lab_list = ['ccle', 'edgren', 'gcsi', 'gray', 'uhn']

    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        your_list = list(reader)
        for n in range(1,5):
            for lab in lab_list:
                new_csv_list = [0, 0, 0]
                for line in your_list:
                    if len(line) > 0:
                        if line[0] == str(n) and line[1] == lab:
                            # do cell average for n
                            new_csv_list[0] += int(line[5])
                            new_csv_list[1] += int(line[6])
                            new_csv_list[2] += int(line[7])

                new_acc = new_csv_list[0] / new_csv_list[1]
                new_prec = new_csv_list[0] / (new_csv_list[0] + new_csv_list[2])
                new_csv_list.insert(0, new_prec)
                new_csv_list.insert(0, new_acc)
                new_csv_list.insert(0, lab)
                new_csv_list.insert(0, n)
                csv_list.append(new_csv_list)

    with open(save_dir + 'cell_line_averages.csv', "w", newline="") as f2:
        writer = csv.writer(f2)
        writer.writerows(csv_list)
    return
