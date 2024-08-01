import pandas as pd
import numpy as np
from scipy.stats import norm, linregress
from typing import List
from collections import namedtuple
from collections import defaultdict


ASYMPTOMATIC_SHARE = 0.35
DOUBLING_PERIOD_D = 3
DEBUG = None


def logit(x):
    return np.log(x / (1 - x))


def logistic(x):
    return 1 / (1 + np.exp(-x))


def get_studies():
    df_op_lu = pd.read_csv(
        "data/2024-06-17-swab-sensitivity/lu_op_ct_mgs.tsv",
        sep="\t",
        skiprows=1,
    )  # Data obtained from Table S1.
    df_op_lu.rename(
        columns={"SCV-2 Relative Abundance": "scv2_ra", "Ct value": "scv2_ct"},
        inplace=True,
    )
    df_op_lu[["patient_status", "swab_type", "Study"]] = [
        "Inpatient",
        "op",
        "Lu et al. 2021",
    ]

    df_np_rodriguez = pd.read_csv(
        "data/2024-06-17-swab-sensitivity/rodriguez_np_ct_mgs.csv", sep=";"
    )  # Data sent to us by authors.
    rodriguez_patient_status_dict = {
        "Hospit": "Inpatient",
        "Out_Patient": "Outpatient",
        "Intensive_Care": "ICU",
    }
    df_np_rodriguez["patient_status"] = df_np_rodriguez["Group"].replace(
        rodriguez_patient_status_dict
    )
    df_np_rodriguez["scv2_ra"] = (
        df_np_rodriguez["Reads_2019_CoV"] / df_np_rodriguez["Reads_Total"]
    )

    df_np_rodriguez.rename(columns={"CoV_Ct_number": "scv2_ct"}, inplace=True)
    df_np_rodriguez[["swab_type", "Study"]] = ["np", "Rodriguez et al. 2021"]

    df_np_babiker = pd.read_csv(
        "data/2024-06-17-swab-sensitivity/babiker_np_ct_mgs.tsv",
        sep="\t",
        skiprows=1,
    )  # Data obtained from table S2
    df_np_babiker.rename(
        columns={
            "SARS-CoV-2 RT-PCR Ct": "scv2_ct",
            "SARS-CoV-2 RA": "scv2_ra",
            "Inpatient/ED vs. Outpatient": "patient_status",
        },
        inplace=True,
    )
    df_np_babiker["scv2_ct"] = (
        df_np_babiker["scv2_ct"].replace(",", ".", regex=True).astype(float)
    )
    df_np_babiker["patient_status"] = df_np_babiker["patient_status"].apply(
        lambda x: x if x in ["Inpatient", "Outpatient"] else "Unknown"
    )
    # The data uses . to represent missing data. Set this column to integers, while at the same time mapping missing data to NA.
    df_np_babiker["days_from_onset"] = (
        df_np_babiker["Day of Testing Relative to Symptom Onset"]
        .replace(".", "-1")
        .astype(int)
        .replace(-1, "NA")
    )
    df_np_babiker[["swab_type", "Study"]] = ["np", "Babiker et al. 2020"]

    df_np_mostafa = pd.read_csv(
        "data/2024-06-17-swab-sensitivity/mostafa-np-ra-ct.tsv", sep="\t"
    )  # Data obtained from Table S2.
    mostafa_severity_dict = {
        1: "Required\nventilator",
        2: "ICU",
        3: "Inpatient",
        4: "Outpatient",
        0: "Unknown",
    }
    df_np_mostafa.rename(
        columns={
            "SARS-CoV-2 RT-PCR Ct value": "scv2_ct",
            "CosmosID Proportion Mapped to SARS-CoV-2": "scv2_ra",
        },
        inplace=True,
    )
    df_np_mostafa["Severity index"] = df_np_mostafa["Severity index"].replace("–", 0)
    df_np_mostafa["patient_status"] = (
        df_np_mostafa["Severity index"].astype(int).replace(mostafa_severity_dict)
    )
    # There is no information of why some patients only have "<7" as their days from onset. We set it to 3.5 (the average of 1-6 days.)
    df_np_mostafa["days_from_onset"] = df_np_mostafa["No. of days from onset"].replace(
        {"–": "NA", "<7": "3.5"}
    )
    # Drop samples unless we have both qPCR and MGS detection
    df_np_mostafa = df_np_mostafa[df_np_mostafa["COVID-19-positive"] == True]
    df_np_mostafa = df_np_mostafa[df_np_mostafa["scv2_ct"] != "–"]
    df_np_mostafa["scv2_ct"] = df_np_mostafa["scv2_ct"].astype(float)
    df_np_mostafa[["swab_type", "Study"]] = ["np", "Mostafa et al. 2020"]

    study_dfs = {
        "Lu et al. 2021": df_op_lu,
        "Babiker et al. 2020": df_np_babiker,
        "Mostafa et al. 2020": df_np_mostafa,
        "Rodriguez et al. 2021": df_np_rodriguez,
    }
    return study_dfs


def get_study_ras():
    studies = get_studies()
    study_ras = {}
    for title, df in studies.items():
        study_ras[title] = df["scv2_ra"].tolist()
    return study_ras


def get_composite_ras():
    study_ras = get_study_ras().values()
    composite_swab_ras = sum(study_ras, [])
    return composite_swab_ras


def get_adjusted_symp_and_asymp_composite_ras():
    study_dfs = get_studies().values()
    composite_df = pd.concat(study_dfs)
    composite_df = composite_df[
        composite_df["patient_status"].isin(["Inpatient", "Outpatient"])
    ]
    zero_ras = composite_df[composite_df["scv2_ra"] == 0]["scv2_ra"].tolist()
    symptom_status_dfs = adjust_cts(composite_df)
    symptom_status_ras = defaultdict(list)
    for symptom_status in ["Asymptomatic", "Symptomatic"]:
        df = symptom_status_dfs[symptom_status]
        df = df[df["scv2_ra"] != 0]
        df = adjust_rel_abun(df)
        ras = df["adjusted_scv2_ra"].tolist() + zero_ras
        symptom_status_ras[symptom_status] = ras

    return symptom_status_ras


def get_asymptomatic_factor():
    # https://doi.org/10.1371/journal.pone.0270694
    long_2020_asymptomatic_delta = 0  # "The initial CT values for 37 asymptomatic individuals and 37 symptomatic patients appeared similar" https://www.nature.com/articles/s41591-020-0965-6
    lee_2020_asymptomatic_delta = 0  # "There were no significant differences in CT values between asymptomatic and symptomatic (including presymptomatic) patients." 10.1001/jamainternmed.2020.3862
    yang_2023_asymptomatic_delta = 0.99  # Extracted from supplement figure 4D https://doi.org/10.1016/S2666-5247(23)00139-8

    hall_asymptomatic_ct_median = 29.9
    hall_symptomatic_ct_median = 21.8
    hall_asymptomatic_delta = hall_asymptomatic_ct_median - hall_symptomatic_ct_median
    ASYMPTOMATIC_ADJUSTMENT_FACTOR = (
        hall_asymptomatic_delta
        + long_2020_asymptomatic_delta
        + lee_2020_asymptomatic_delta
        + yang_2023_asymptomatic_delta
    ) / 4

    return ASYMPTOMATIC_ADJUSTMENT_FACTOR


def adjust_cts(df):
    ASYMPTOMATIC_ADJUSTMENT_FACTOR = get_asymptomatic_factor()

    np_data = pd.read_csv(
        "data/2024-06-17-swab-sensitivity/2024-06-18-np-nasal-ct.tsv",
        sep="\t",
        skiprows=1,
    )
    np_means = np_data.mean()

    NP_ADJUSTMENT_FACTOR = np_means.mean()
    goodall_data = pd.read_csv(
        "data/2024-06-17-swab-sensitivity/goodall-op-nasal-ct.tsv",
        sep="\t",
        skiprows=2,
        header=None,
    )
    OP_ADJUSTMENT_FACTOR = goodall_data[0].mean()
    if DEBUG:
        print(f"Goodall mean (OP): {op_goodall_mean}")
        print(f"NP mean: {NP_ADJUSTMENT_FACTOR}")

    df["adjusted_scv2_ct"] = df["scv2_ct"]
    # Subtract the adjustment factors from the CT values (NP_ADJUSTMENT_FACTOR is negative, so it increases the CT values)
    df.loc[df["swab_type"] == "np", "adjusted_scv2_ct"] -= NP_ADJUSTMENT_FACTOR
    df.loc[df["swab_type"] == "op", "adjusted_scv2_ct"] -= OP_ADJUSTMENT_FACTOR
    df_symptomatic = df.copy()
    df_asymptomatic = df.copy()
    df_asymptomatic["adjusted_scv2_ct"] += ASYMPTOMATIC_ADJUSTMENT_FACTOR

    symptom_status_dfs = {
        "Asymptomatic": df_asymptomatic,
        "Symptomatic": df_symptomatic,
    }

    return symptom_status_dfs


def adjust_rel_abun(composite_df):
    composite_df = composite_df.copy()
    composite_df.loc[:, "scv2_ra_logit"] = composite_df["scv2_ra"].apply(logit)

    slope, intercept, r_value, p_value, std_err = linregress(
        composite_df["scv2_ct"], composite_df["scv2_ra_logit"]
    )
    composite_df["adjusted_scv2_ra_logit"] = (
        intercept + slope * composite_df["adjusted_scv2_ct"]
    )
    residuals = composite_df["scv2_ra_logit"] - (
        intercept + slope * composite_df["scv2_ct"]
    )

    sigma_squared = np.var(residuals, ddof=2)

    sigma = np.sqrt(sigma_squared)

    noise = np.random.normal(loc=0, scale=sigma, size=len(composite_df))

    composite_df["adjusted_scv2_ra_logit_with_noise"] = (
        composite_df["adjusted_scv2_ra_logit"] + noise
    )
    composite_df["adjusted_scv2_ra"] = composite_df[
        "adjusted_scv2_ra_logit_with_noise"
    ].apply(logistic)
    return composite_df


def get_logit_normal_samples(ras):
    ra_values = np.array(ras)
    zero_share = (ra_values == 0).mean()
    ra_values = ra_values[ra_values != 0]
    logit_ra_values = logit(ra_values)
    mean, std = np.mean(logit_ra_values), np.std(logit_ra_values)
    norm_dist = norm(loc=mean, scale=std)
    logit_samples = norm_dist.rvs(size=int(100000 * (1 - zero_share)))
    samples = logistic(logit_samples)
    samples = np.append(samples, np.zeros(int(100000 * zero_share)))
    np.random.shuffle(samples)
    return samples


def get_adjusted_composite_ras():
    symptom_status_ras = get_adjusted_symp_and_asymp_composite_ras()

    asymptomatic_samples = get_logit_normal_samples(symptom_status_ras["Asymptomatic"])
    symptomatic_samples = get_logit_normal_samples(symptom_status_ras["Symptomatic"])

    adjusted_composite_ras = np.random.choice(
        asymptomatic_samples, round(len(asymptomatic_samples) * ASYMPTOMATIC_SHARE)
    )

    adjusted_composite_ras = np.append(
        adjusted_composite_ras,
        np.random.choice(
            symptomatic_samples,
            round(len(symptomatic_samples) * (1 - ASYMPTOMATIC_SHARE)),
        ),
    )

    return adjusted_composite_ras


adjusted_composite_ras = get_adjusted_composite_ras()

# Save adjusted composite RAs as a TSV file
output_df = pd.DataFrame({"relative_abundance": adjusted_composite_ras})
output_df.to_csv("data/adjusted_composite_ras.tsv", sep="\t", index=False)
