#! /usr/bin/env python3

import subprocess
import re
import csv
import matplotlib.pyplot as plt
import numpy as np

EXE_NAME = "usi_itp_ocl"
USI_NAME = "image.dat"

cpu_time_pat = re.compile(r'(\w+):.*CPU\D*([1-9]\d*\.\d*|0\.\d*)')
ocl_time_pat = re.compile(r'(\w+ocl)\D*([1-9]\d*\.\d*|0\.\d*)')

cpu_time = {}
ocl_time = {}
output = subprocess.check_output([EXE_NAME, USI_NAME], stderr=subprocess.STDOUT).decode()
lines = output.splitlines()
for line in lines:
    m = cpu_time_pat.match(line)
    if(m):
        cpu_time[m.group(1)] = float(m.group(2))
    
    m = ocl_time_pat.match(line)
    if(m):
        ocl_time[m.group(1)] = float(m.group(2))

print(cpu_time)
print(ocl_time)

acc_ratio = {}

COL_NAMES = ["method", "CPU time (ms)", "OpenCL GPU time (ms)", "accelerate ratio"]
with open('benchmark.csv','w', newline='') as f:
    f_csv = csv.writer(f)
    f_csv.writerow(COL_NAMES)

    for method_cpu in cpu_time.keys():
        for method_ocl in ocl_time.keys():
            if method_ocl == method_cpu + '_ocl':
                acc_ratio[method_cpu] = cpu_time[method_cpu] / ocl_time[method_ocl]
                f_csv.writerow([method_cpu,
                                '%.2f' % cpu_time[method_cpu], 
                                '%.2f' % ocl_time[method_ocl], 
                                '%.2f' % acc_ratio[method_cpu]])
print(acc_ratio)

plt.figure('evaluation', figsize=(9, 7))

measures = {'CPU': cpu_time, 'OpenCL GPU': ocl_time, 'Accelerate Ratio': acc_ratio}
for i, (name, data) in enumerate(measures.items()):
    ax = plt.subplot2grid((3, 1), (i, 0))
    ax.set_title(name)
    color = 'r' if name == 'Accelerate Ratio' else None
    ax.barh(list(data.keys()), list(data.values()), color=color, alpha=0.6)
    # xpos = np.arange(len(data))
    # ax.barh(xpos, list(data.values()), color=color, alpha=0.6)
    # plt.yticks(xpos, list(data.keys()))
    for k,v in data.items():
        ax.text(v, k, '%.2f' % v, ha='left', va= 'center')

plt.tight_layout()
plt.show()