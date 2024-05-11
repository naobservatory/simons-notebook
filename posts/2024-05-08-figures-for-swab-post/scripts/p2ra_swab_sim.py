#!/usr/bin/env python3
import pandas as pd
import numpy as np
from collections import defaultdict
import random
import math
import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt

def simulate_p2ra_new_many(samples, sample_populations=[100],  n_simulations=1000) -> pd.DataFrame:
    results = defaultdict(list)
    for sample_pop in sample_populations:
        for _ in range(n_simulations):
            results[sample_pop].append(simulate_p2ra_new( samples, sample_pop))
    for key, values in results.items():
        results[key] = sorted(values)
    df = pd.DataFrame(results)
    return df


def simulate_p2ra_new(samples, sample_pop):
    target_incidence_p_w = 0.01
    shedding_duration_w = 2

    prevalence = target_incidence_p_w * shedding_duration_w
    n_sick = np.random.poisson(sample_pop * prevalence)
    print(n_sick)
    if n_sick == 0:
        return 0
    ra_sick = 0
    for _ in range(n_sick):
        observation = np.random.choice(samples)
        ra_sick += observation
    ra_sick = ra_sick / n_sick
    relative_abundance = n_sick / sample_pop * ra_sick
    #if relative_abundance > 1:
    #    print(n_sick, ra_sick, relative_abundance)
    #    print("Large relative abundance.")
    return relative_abundance



def get_p2ra_dfs():
    df_np_babiker, df_np_rodriguez, df_op_lu, df_np_mostafa = return_studies()
    df_wastewater = pd.read_csv('../data/2024-05-07-fits.tsv', sep='\t')
    df_op_lu_ras = df_op_lu['scv2_ra'].dropna().tolist()
    df_np_babiker_ras = df_np_babiker['scv2_ra'].dropna().tolist()
    df_np_mostafa_ras = df_np_mostafa['scv2_ra'].dropna().tolist()
    df_np_rodriguez_ras = df_np_rodriguez['scv2_ra'].dropna().tolist()
    combined_swab_ras = df_op_lu_ras + df_np_babiker_ras + df_np_mostafa_ras + df_np_rodriguez_ras
    df_rothman_ras = df_wastewater[(df_wastewater['study'] == 'rothman') & (df_wastewater['location'] == 'Overall') & (df_wastewater['pathogen'] == 'sars_cov_2')]['ra_at_1in100'].tolist()
    df_crits_christoph_ras = df_wastewater[(df_wastewater['study'] == 'crits_christoph') & (df_wastewater['location'] == 'Overall') & (df_wastewater['pathogen'] == 'sars_cov_2')]['ra_at_1in100'].tolist()



    swab_studies = ["Lu et al. 2021", "Babiker et al. 2020", "Mostafa et al. 2020", "Rodriguez et al. 2021", "Combined"]
    n_swabs = [50, 100, 200]
    swab_df = pd.DataFrame()

    def logit(x):
        return np.log(x / (1 - x))

    def logistic(x):
        return 1 / (1 + np.exp(-x))

    raw_distributions = []

    for study, positive_ras in zip(swab_studies, [df_op_lu_ras, df_np_babiker_ras, df_np_mostafa_ras, df_np_rodriguez_ras, combined_swab_ras]):
        ra_values = np.array(positive_ras)
        ra_values = logit(ra_values)

        mean, std = np.mean(ra_values), np.std(ra_values)
        norm_dist = scipy.stats.norm(loc=mean, scale=std)
        logit_samples = norm_dist.rvs(size=100000)
        samples = logistic(logit_samples)
        raw_distributions.append(samples)

        df = simulate_p2ra_new_many(samples, n_swabs, n_simulations=10000)

        df["Study"] = study

        swab_df = pd.concat([swab_df, df.melt(id_vars=["Study"], value_vars=n_swabs, var_name="Sample Size",
        value_name="Relative Abundance")])

    fig, ax = plt.subplots(1, 1, figsize=(8, 6), dpi=600)
    #for distribution in raw_distributions:
    #    log_distribution = np.log10(distribution)
    #    sns.kdeplot(log_distribution, ax=ax)
    plt.tight_layout()
    plt.savefig('../fig/test-plot.png', dpi=600)
    plt.clf()

    swab_df.reset_index(drop=True, inplace=True)
    df_ww = pd.DataFrame({
        'Relative Abundance': df_rothman_ras + df_crits_christoph_ras,
        'Study': ['Rothman et al. 2021'] * len(df_rothman_ras) + ['Crits-Christoph et al. 2021'] * len(df_crits_christoph_ras)
    })

    return swab_df, df_ww

def return_fig():
    plot_df_swab, plot_df_ww = get_p2ra_dfs()
    plot_df_ww['Relative Abundance'] = np.log10(plot_df_ww['Relative Abundance'])

    plot_df_swab['Relative Abundance'] = np.log10(plot_df_swab['Relative Abundance'])

    fig, axs = plt.subplots(2, 1, figsize=(8, 6), dpi=600, gridspec_kw={'height_ratios': [3, 1]})

    fig.subplots_adjust(hspace=0.4)
    colors = sns.color_palette("viridis", 3)

    sns.violinplot(x='Relative Abundance', y='Study', ax=axs[0], hue='Sample Size', palette=colors, data=plot_df_swab, inner=None, linewidth=0, bw_adjust=1)
    sns.violinplot(x='Relative Abundance', y='Study', ax=axs[1], color='#aedc31', data=plot_df_ww, inner=None, linewidth=0, bw_adjust=1, width=0.5)

    axs[0].set_title('A', x=-0.29, y=0.9)
    axs[1].set_title('B', x=-0.29, y=0.9)

    for ax in axs:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_ylabel('')

        y_labels = [label.get_text().rsplit(' ', 1)[0] + '\n' + label.get_text().rsplit(' ', 1)[1] if ' ' in label.get_text() else label.get_text() for label in ax.get_yticklabels()]
        ax.set_yticklabels(y_labels)

        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=True)


    xmin, xmax = axs[0].get_xlim()

    for x in np.arange(math.ceil(xmin), 1, 1):
        axs[0].axvline(x=x, color='black', linestyle='--', linewidth=0.3, alpha=0.2, zorder=-1)
        axs[1].axvline(x=x, color='black', linestyle='--', linewidth=0.3, alpha=0.2, zorder=-1)

    current_xticks = axs[1].get_xticks()
    for ax in axs:
        ax.set_xlim(xmin, 1)
        ax.set_xticklabels(['$10^{{{}}}$'.format(int(x)) if x != 0 else '1' for x in current_xticks])


    axs[0].legend(title='Number of swabs', loc='center left', bbox_to_anchor=(0.9, 0.5), ncol=1, labels=['50', '100', '200'], frameon=False)
    plt.tight_layout()
    plt.savefig('../fig/swab-ww-ra_1.png', dpi=600)


def return_studies():
    df_op_lu = pd.read_csv('../data/lu_throat_ct_mgs.tsv', sep='\t', skiprows=1)
    df_op_lu.rename(columns={'SCV-2 Relative Abundance': 'scv2_ra', 'Ct value': 'scv2_ct'}, inplace=True)
    df_op_lu['patient_status'] = 'Inpatient'  # All inpatients for Lu et al.
    df_op_lu["swab_type"] = "op"

    df_np_babiker = pd.read_csv('../data/babiker_np_ct_mgs.tsv', sep='\t', skiprows=1)
    df_np_babiker.rename(columns={'SARS-CoV-2 RT-PCR Ct': 'scv2_ct', 'SARS-CoV-2 RA': 'scv2_ra', 'Inpatient/ED vs. Outpatient': 'patient_status'}, inplace=True)
    df_np_babiker['scv2_ct'] = df_np_babiker['scv2_ct'].replace(',', '.', regex=True).astype(float)
    df_np_babiker['patient_status'] = df_np_babiker['patient_status'].apply(lambda x: x if x in ['Inpatient', 'Outpatient'] else 'Unknown')
    df_np_babiker["swab_type"] = "np"

    df_np_mostafa = pd.read_csv('../data/mostafa_np_scv2_ct_mgs.tsv', sep='\t', skiprows=1)
    mostafa_severity_dict = {
        1: "Required\nventilator",
        2: "ICU",
        3: "Inpatient",
        4: "Outpatient",
        0: "Unknown"
    }
    df_np_mostafa.rename(columns={'SARS-CoV-2 RT-PCR Ct': 'scv2_ct', 'SARS-CoV-2 RA': 'scv2_ra'}, inplace=True)
    df_np_mostafa['patient_status'] = df_np_mostafa['Severity Index'].astype(int).replace(mostafa_severity_dict)
    df_np_mostafa["swab_type"] = "np"


    df_np_rodriguez = pd.read_csv('../data/rodriguez_np_ct_mgs.csv', sep=';')
    rodriguez_patient_status_dict = {
        "Hospit": "Inpatient",
        "Out_Patient": "Outpatient",
        "Intensive_Care": "ICU"
    }

    df_np_rodriguez.rename(columns={"CoV_Ct_number": "scv2_ct"}, inplace=True)
    df_np_rodriguez['patient_status'] = df_np_rodriguez['Group'].replace(rodriguez_patient_status_dict)
    df_np_rodriguez["scv2_ra"] = df_np_rodriguez["Reads_2019_CoV"] / df_np_rodriguez["Reads_Post_trimming"]
    df_np_rodriguez = df_np_rodriguez[df_np_rodriguez["scv2_ra"] != 0] #
    df_np_rodriguez["swab_type"] = "np"

    return [df_np_babiker, df_np_rodriguez, df_op_lu, df_np_mostafa]

def rpkm_to_ra(df):
    genome_lengths_in_kb = {
        "Influenza A": 13.5,
        "Influenza B": 14.5,
        "Metapneumovirus": 13.3,
        "Rhinovirus": 7.2,
        "Parainfluenzavirus 1": 15.5,
        "Parainfluenzavirus 3": 15.5,
        "Parainfluenzavirus 4": 17.4,
        "Respiratory Syncytial Virus": 15.2
    }


    for virus, virus_length in genome_lengths_in_kb.items():
        df.loc[df["virus"] == virus, "relative_abundance"] = df.loc[df["virus"] == virus, "RPKM"] * virus_length / 1e6

    return df


if __name__ == "__main__":
    return_fig()

