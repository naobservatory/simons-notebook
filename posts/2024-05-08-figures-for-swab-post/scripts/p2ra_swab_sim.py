#!/usr/bin/env python3
import pandas as pd
import numpy as np
from collections import defaultdict
import random
import math
import seaborn as sns
import matplotlib.pyplot as plt

def simulate_p2ra_new_many(shedding_values=[0.01], sample_populations=[100],  n_simulations=1000) -> pd.DataFrame:
    results = defaultdict(list)
    for sample_pop in sample_populations:
        for _ in range(n_simulations):
            results[sample_pop].append(simulate_p2ra_new(sample_pop, shedding_values))
    for key, values in results.items():
        results[key] = sorted(values)
    df = pd.DataFrame(results)
    return df


def simulate_p2ra_new(sample_pop=100, shedding_values=[0.01]):
    target_incidence_p_w = 0.01
    shedding_duration_w = 2
    sigma_shedding_values = 0.05
    epsilon = 0.000001
    prevalence = target_incidence_p_w * shedding_duration_w
    n_sick = np.random.poisson(sample_pop * prevalence)
    ra_sicks = get_inputs_biased(shedding_values, sigma_shedding_values, epsilon)

    ra_sick = 0
    for _ in range(n_sick):
        ra_sick += ra_sicks[random.randint(0, len(ra_sicks) - 1)]

    relative_abundance = n_sick / sample_pop * ra_sick
    return relative_abundance


def get_inputs_biased(shedding_values, sigma_shedding_values, epsilon):
    empirical_values = shedding_values
    sigma = sigma_shedding_values
    if sigma < epsilon:
        return empirical_values

    # Don't want to bias each mean independently, bias all of them
    # together with log-normally distributed noise. The geometric mean
    # (and median) of the noise is zero, and the standard deviation is
    # provided by the user.

    bias = np.random.lognormal(0, sigma)

    adjusted_values = [empirical_value * bias for empirical_value in empirical_values]
    return adjusted_values



def format_func(value, tick_number):
    return r'$10^{{{}}}$'.format(int(value))

def return_studies():
    df_op_lu = pd.read_csv('../data/lu_throat_ct_mgs.tsv', sep='\t', skiprows=1)
    df_op_lu.rename(columns={'SCV-2 Relative Abundance': 'scv2_ra', 'Ct value': 'scv2_ct'}, inplace=True)
    df_op_lu['patient_status'] = 'Inpatient'
    df_np_babiker = pd.read_csv('../data/babiker_np_ct_mgs.tsv', sep='\t', skiprows=1)
    df_np_babiker.rename(columns={'SARS-CoV-2 RT-PCR Ct': 'scv2_ct', 'SARS-CoV-2 RA': 'scv2_ra', 'Inpatient/ED vs. Outpatient': 'patient_status'}, inplace=True)
    df_np_mostafa = pd.read_csv('../data/mostafa_np_scv2_ct_mgs.tsv', sep='\t', skiprows=1)
    df_np_rodriguez = pd.read_csv('../data/rodriguez_np_ct_mgs.csv', sep=';')

    mostafa_severity_dict = {
        1: "Required\nventilator",
        2: "ICU",
        3: "Inpatient",
        4: "Outpatient",
        0: "Unknown"
    }

    df_np_mostafa.rename(columns={'SARS-CoV-2 RT-PCR Ct': 'scv2_ct', 'SARS-CoV-2 RA': 'scv2_ra'}, inplace=True)
    df_np_mostafa['patient_status'] = df_np_mostafa['Severity Index'].astype(int).replace(mostafa_severity_dict)

    rodriguez_patient_status_dict = {
        "Hospit": "Inpatient",
        "Out_Patient": "Outpatient",
        "Intensive_Care": "ICU"
    }

    df_np_rodriguez['patient_status'] = df_np_rodriguez['Group'].replace(rodriguez_patient_status_dict)
    df_np_rodriguez["scv2_ra"] = df_np_rodriguez["Reads_2019_CoV"] / df_np_rodriguez["Reads_Post_trimming"]
    df_np_rodriguez = df_np_rodriguez[df_np_rodriguez["scv2_ra"] != 0] # Dropping Patient_066 because 0 ra breaks gmean

    df_np_rodriguez.rename(columns={"CoV_Ct_number": "scv2_ct"}, inplace=True)

    df_np_babiker['scv2_ct'] = df_np_babiker['scv2_ct'].replace(',', '.', regex=True).astype(float)

    df_op_lu['patient_status'] = 'Inpatient'  # All inpatients for Lu et al.
    df_np_babiker['patient_status'] = df_np_babiker['patient_status'].apply(lambda x: x if x in ['Inpatient', 'Outpatient'] else 'Unknown')

    df_op_lu["swab_type"] = "op"
    df_np_babiker["swab_type"] = "np"
    df_np_rodriguez["swab_type"] = "np"
    df_np_mostafa["swab_type"] = "np"
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
    for study, positive_ras in zip(swab_studies, [df_op_lu_ras, df_np_babiker_ras, df_np_mostafa_ras, df_np_rodriguez_ras, combined_swab_ras]):
        df = simulate_p2ra_new_many(positive_ras, n_swabs, n_simulations=2000)

        df["Study"] = study

        swab_df = pd.concat([swab_df, df.melt(id_vars=["Study"], value_vars=n_swabs, var_name="Sample Size",
        value_name="Relative Abundance")])
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

    sns.violinplot(x='Relative Abundance', y='Study', ax=axs[0], hue='Sample Size', palette=colors, data=plot_df_swab, inner=None, linewidth=0, bw_adjust=0.3)
    sns.violinplot(x='Relative Abundance', y='Study', ax=axs[1], color='#aedc31', data=plot_df_ww, inner=None, linewidth=0, bw_adjust=0.3, width=0.5)


    #for study in swab_studies:
    #    for sample_size in n_swabs:
    #        zero_share = len(swab_df[(swab_df['Study'] == study) & (swab_df['Sample Size'] == sample_size) & (swab_df['Relative Abundance'] == 0)]) / len(swab_df[(swab_df['Study'] == study) & (swab_df['Sample Size'] == sample_size)])
    #        print(study, sample_size, zero_share)


    # Adding titles to the subplots
    axs[0].set_title('A', x=-0.29, y=0.9)
    axs[1].set_title('B', x=-0.29, y=0.9)

    for ax in axs:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_ylabel('')

        # Adjusting y labels to have a newline before the year
        y_labels = [label.get_text().rsplit(' ', 1)[0] + '\n' + label.get_text().rsplit(' ', 1)[1] if ' ' in label.get_text() else label.get_text() for label in ax.get_yticklabels()]
        ax.set_yticklabels(y_labels)

        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=True)

    # Retrieve current x-ticks from the first axis to apply to both
    xmin, xmax = axs[0].get_xlim()
    print(xmin, xmax)

    for x in np.arange(math.ceil(xmin), 1, 1):
        axs[0].axvline(x=x, color='black', linestyle='--', linewidth=0.3, alpha=0.2, zorder=-1)
        axs[1].axvline(x=x, color='black', linestyle='--', linewidth=0.3, alpha=0.2, zorder=-1)

    current_xticks = axs[1].get_xticks()
    for ax in axs:
        ax.set_xlim(xmin, 1)
        ax.set_xticklabels(['$10^{{{}}}$'.format(int(x)) if x != 0 else '1' for x in current_xticks])

    print(current_xticks)
    axs[0].legend(title='Number of swabs', loc='center left', bbox_to_anchor=(0.9, 0.5), ncol=1, labels=['50', '100', '200'], frameon=False)
    plt.tight_layout()
    plt.savefig('../fig/swab-ww-ra_1.png', dpi=600)

    plt.show()

if __name__ == "__main__":
    return_fig()

