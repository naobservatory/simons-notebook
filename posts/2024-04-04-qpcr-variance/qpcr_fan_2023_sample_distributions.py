# Importing packages
# Does not work.
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

colors = sns.color_palette("tab10")


def sample_distributions(list_of_means, list_of_sds):
    aggregate_means_sds = []
    for means, sds in zip(list_of_means, list_of_sds):
        N_PICKS = 500
        samples = []
        for mean, sd in zip(means, sds):
            try:
                samples.append(np.random.normal(mean, sd, N_PICKS))
            except:
                print(mean, sd)
        aggregate_means_sds.append([np.mean(samples), np.std(samples)])
    return aggregate_means_sds



with open("data/fan_2023_supplement_Table_S1.tsv") as f:
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    qpcr_target_na_kit_machine = [
        ["orf1ab", "different", "different"],
        ["N", "different", "different"],
        ["orf1ab", "Tianlong", "different"],
        ["N", "Tianlong", "different"],
        ["orf1ab", "Tianlong", "ABI7500"],
        ["N", "Tianlong", "ABI7500"],
        ["orf1ab", "DaAn", "different"],
        ["N", "DaAn", "different"]
    ]

    lines = f.readlines()
    lines = lines[2:] # skipping first two lines
    table_iterator = -1
    for line in lines:
        bits = line.split("\t")

        if bits == ["\n"]:
            continue


        if "different" in bits[0]:

            if table_iterator >= 0: # returning values from previous subtabpe
                means = [means_2_e3, means_1_e3, means_5_e2, means_2_e2]
                sds = [sds_2_e3, sds_1_e3, sds_5_e2, sds_2_e2]
                aggregate_means_sds = sample_distributions(means, sds)
                print("hello", aggregate_means_sds)

            table_iterator += 1  # new subtable, resetting.
            qpcr_info = qpcr_target_na_kit_machine[table_iterator]
            color_iterator = 0
            increment = -0.3
            means_2_e3, means_1_e3, means_5_e2, means_2_e2 = [], [], [], []
            sds_2_e3, sds_1_e3, sds_5_e2, sds_2_e2 = [], [], [], []
            continue



        increment += 0.1
        qpcr_kit, mean_2e3, sd_2e3, n_labs_2e3, _, mean_1e3, sd_1e3, n_labs_1e3, _, mean_5e2, sd_5e2, n_labs_5e2, _, mean_2e2, sd_2e2,n_labs_2e2 = bits
        means = [float(mean_2e3), float(mean_1e3), float(mean_5e2), float(mean_2e2)]
        sds = [float(sd_2e3), float(sd_1e3), float(sd_5e2), float(sd_2e2)]

        means_2_e3.append(mean_2e3)
        means_1_e3.append(mean_1e3)
        means_5_e2.append(mean_5e2)
        means_2_e2.append(mean_2e2)

        sds_2_e3.append(sd_2e3)
        sds_1_e3.append(sd_1e3)
        sds_5_e2.append(sd_5e2)
        sds_2_e2.append(sd_2e2)





        continue
