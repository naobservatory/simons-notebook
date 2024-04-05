
# Importing packages
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

colors = sns.color_palette("tab10")


def plot_dot_and_sd(ax, color_iterator, increment, y_means, y_sds,qpcr_kit):
    x_values = [x + increment for x in range(len(y_means))]
    ax.errorbar(x_values, y_means, y_sds,  color=colors[color_iterator], alpha=0.5, label=qpcr_kit, fmt='-o', linestyle='None')

    # set x axis values to concentrations
    ax.set_xticks(range(len(y_means)))
    ax.set_xticklabels(["2e3", "1e3", "5e2", "2e2"])


with open("data/fan_2023_supplement_Table_S1.tsv") as f:
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))


    lines = f.readlines()
    target_sub_table = False
    table_iterator = -1

    for line in lines:
        bits = line.split("\t")
        if bits == ["\n"]:
            continue

        if bits[0].count("different") == 1:
            target_sub_table = True
            table_iterator += 1
            color_iterator = 0
            increment = -0.3
            continue
        if bits[0].count("different") == 2:
            target_sub_table = False
            continue
        if not target_sub_table:
            continue

        if target_sub_table:
            color_iterator +=1
            increment += 0.1
            qpcr_kit, mean_2e3, sd_2e3, n_labs_2e3, _, mean_1e3, sd_1e3, n_labs_1e3, _, mean_5e2, sd_5e2, n_labs_5e2, _, mean_2e2, sd_2e2, n_labs_2e2 = bits
            qpcr_kit = qpcr_kit.strip()

            # Replace empty strings with a default value before converting to float
            means = [float(mean) if mean else 0.0 for mean in [mean_2e3, mean_1e3, mean_5e2, mean_2e2]]
            sds = [float(sd) if sd else 0.0 for sd in [sd_2e3, sd_1e3, sd_5e2, sd_2e2]]

            plot_dot_and_sd(axs[table_iterator], color_iterator, increment, means, sds, qpcr_kit)

    # title for ax0 = ORF1ab
    # title for ax1 = N gene
    axs[0].set_title("ORF1ab gene")
    axs[1].set_title("N gene")
    plt.suptitle("Ct values for different RT-PCR kits, with NA kit and PCR machine being the same.")
    axs[1].legend(loc="lower right", title="RT-PCR kit")
    # add horitzontal lines
    #for i in range()
    #    ax.axhline(30, color="black", linestyle="--")
    plt.tight_layout()
    for i in range (28,40,2):
        axs[0].axhline(i, color="grey", linestyle="--", linewidth=0.5)
        axs[1].axhline(i, color="grey", linestyle="--")


    plt.show()
    plt.clf()


