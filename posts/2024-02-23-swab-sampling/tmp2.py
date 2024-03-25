#!/usr/bin/env python3

import numpy as np
import pandas as pd
import math
import random
import matplotlib.pyplot as plt
from collections import defaultdict
def simulate_detection(positive_ras):
    FRAC_SICK = 0.01  
    SIMULATIONS = 10000
    
    def probabilistic_round(x):
        return int(math.floor(x + random.random()))
    
    def simulate_once(n_swabs):
        n_sick = probabilistic_round(FRAC_SICK*n_swabs)
        if not n_sick:
            return 0
        ra = 0
        for _ in range(n_sick):
            ra = random.choice(positive_ras) 
        return ra / n_swabs
    
    def simulate_many(n_simulations, n_swabs):
        RAs = []
        for _ in range(n_simulations):
            RAs.append(simulate_once(n_swabs))
        RAs.sort()
    
        return [RAs[probabilistic_round(len(RAs)*percentile/100)]
                for percentile in range(100)]
    data = defaultdict(list)
    n_swab_range = [50, 100, 200, 400, 800]
    for n_swabs in n_swab_range: 
        for percentile, ra in enumerate(simulate_many(SIMULATIONS, n_swabs)):
            data[percentile].append(ra)
    df = pd.DataFrame(data).T
    df.columns = n_swab_range 

    return df
df_throat = pd.read_csv('data/lu_throat_ct_mgs.tsv', sep='\t', skiprows=1)
df_throat.rename(columns={'SARS-CoV-2 RT-PCR Ct': 'scv2_ct', 'SARS-CoV-2 RA': 'scv2_ra', 'Inpatient/ED vs. Outpatient': 'patient_status'}, inplace=True)
df_throat['scv2_ct'] = df_throat['scv2_ct'].replace(',', '.', regex=True).astype(float)

df_nasopharyngeal = pd.read_csv('data/babiker_np_ct_mgs.tsv', sep='\t', skiprows=1)
df_nasopharyngeal.rename(columns={'SCV-2 Relative Abundance': 'scv2_ra', 'Ct value': 'scv2_ct'}, inplace=True)

df_throat_ras = df_throat['scv2_ra'].dropna().tolist()
print(sorted(df_throat_ras))
df_nasopharyngeal_ras = df_nasopharyngeal['scv2_ra'].dropna().tolist()
print(sorted(df_nasopharyngeal_ras))
fig, axs = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
n_swab_range = [50, 100, 200, 400, 800]
for i, positive_ras in enumerate([df_throat_ras, df_nasopharyngeal_ras]):
    df= simulate_detection(positive_ras)
    for n_swabs in n_swab_range:
        axs[i].plot(df.index, df[n_swabs], label=f'{n_swabs} swabs')
    axs[i].set_xlabel('Percentile')
    axs[i].set_ylabel('RA')
    axs[i].set_yscale('log')
    # turn on y tick labels for both plots
    axs[i].yaxis.set_tick_params(labelleft=True)
    for x in range (0, 100, 20):
        axs[i].axvline(x=x, color='black', linestyle='--', linewidth=0.5, alpha=0.5)

    for y in range (-1, -10, -1):
        log_y = 10**y
        print(log_y)
        axs[i].axhline(y=log_y, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
    if i == 0:
        axs[i].set_title('Throat')    
    else:
        axs[i].set_title('Nasopharyngeal')


    if i == 1:
        axs[i].legend()


# log y axis
plt.show()
