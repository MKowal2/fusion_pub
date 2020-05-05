import sys
import os
import xlrd
import xlwt
import csv
from glob import glob

path_to_workbook = './tool_outputs_WOC/validated_fusions.xlsx'
workbook_gt = xlrd.open_workbook(path_to_workbook)

# create fusion groundtruth lists
bt474_gt_fusions, bt474_gt_first, bt474_gt_second = return_gt_fusions(workbook_gt,'BT474')
kpl4_gt_fusions, kpl4_gt_first, kpl4_gt_second = return_gt_fusions(workbook_gt,'KPL4')
mcf7_gt_fusions, mcf7_gt_first, mcf7_gt_second = return_gt_fusions(workbook_gt,'MCF7')
skbr3_gt_fusions, skbr3_gt_first, skbr3_gt_second = return_gt_fusions(workbook_gt,'SKBR3')

# save files in columns -> | both genes | first gene | second gene |

bt474 = zip(bt474_gt_fusions,bt474_gt_first,bt474_gt_second)
with open('validated_fusions/bt474_gt.csv', "w") as f:
    writer = csv.writer(f)
    for row in bt474:
        writer.writerow(row)

kpl4 = zip(kpl4_gt_fusions,kpl4_gt_first,kpl4_gt_second)
with open('validated_fusions/kpl4_gt.csv', "w") as f:
    writer = csv.writer(f)
    for row in kpl4:
        writer.writerow(row)

mcf7 = zip(mcf7_gt_fusions,mcf7_gt_first,mcf7_gt_second)
with open('validated_fusions/mcf7_gt.csv', "w") as f:
    writer = csv.writer(f)
    for row in mcf7:
        writer.writerow(row)

skbr3 = zip(skbr3_gt_fusions,skbr3_gt_first,skbr3_gt_second)
with open('validated_fusions/skbr3_gt.csv', "w") as f:
    writer = csv.writer(f)
    for row in skbr3:
        writer.writerow(row)
