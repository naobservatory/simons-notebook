
# Importing packages
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

colors = sns.color_palette("tab10")


def plot_line_and_sd(ax, color_iterator, y_means, y_sds,qpcr_kit):
    lower_sd = [mean - sd for mean, sd in zip(y_means, y_sds)]
    upper_sd = [mean + sd for mean, sd in zip(y_means, y_sds)]

    x_values = range(len(y_means))
    ax.plot(x_values, lower_sd, linestyle='--', color=colors[color_iterator], alpha=0.1)

    ax.plot(x_values, upper_sd, linestyle='--', color=colors[color_iterator], alpha=0.1)
    ax.fill_between(x_values, lower_sd, upper_sd, color=colors[color_iterator], alpha=0.01)
    ax.plot(x_values, y_means, color=colors[color_iterator], label=qpcr_kit)
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
            continue
        if bits[0].count("different") == 2:
            target_sub_table = False
            continue
        if not target_sub_table:
            continue

        if target_sub_table:
            color_iterator +=1
            qpcr_kit, mean_2e3, sd_2e3, n_labs_2e3, _, mean_1e3, sd_1e3, n_labs_1e3, _, mean_5e2, sd_5e2, n_labs_5e2, _, mean_2e2, sd_2e2,n_labs_2e2 = bits
            qpcr_kit = qpcr_kit.strip()

            means = [float(mean_2e3), float(mean_1e3), float(mean_5e2), float(mean_2e2)]
            sds = [float(sd_2e3), float(sd_1e3), float(sd_5e2), float(sd_2e2)]

            plot_line_and_sd(axs[table_iterator], color_iterator, means, sds,qpcr_kit)
            continue

    # title for ax0 = ORF1ab
    # title for ax1 = N gene
    axs[0].set_title("ORF1ab gene")
    axs[1].set_title("N gene")
    plt.suptitle("Ct values for different RT-PCR kits, with NA kit and PCR machine being the same.")
    axs[1].legend(loc="lower right", title="RT-PCR kit")
    plt.tight_layout()

    plt.show()
    plt.clf()



with open("data/fan_2023_supplement_Table_S1.tsv") as f:
    lines = f.readlines()
    target_sub_table = False
    table_iterator = -1
    qpcr_values = []
    for line in lines:
        bits = line.split("\t")
        if bits == ["\n"]:
            continue

        if bits[0].count("different") == 1:
            target_sub_table = True
            table_iterator += 1
            color_iterator = 0
            gene_target = bits[0].split("Ct values for ")[1].split(" ")[0]
            continue
        if bits[0].count("different") == 2:
            target_sub_table = False
            continue
        if not target_sub_table:
            continue




        if target_sub_table:

            qpcr_kit, mean_2e3, sd_2e3, n_labs_2e3, _, mean_1e3, sd_1e3, n_labs_1e3, _, mean_5e2, sd_5e2, n_labs_5e2, _, mean_2e2, sd_2e2,n_labs_2e2 = bits
            mean_2e3 = float(mean_2e3)
            mean_1e3 = float(mean_1e3)
            mean_5e2 = float(mean_5e2)
            mean_2e2 = float(mean_2e2)

            sd_2e3 = float(sd_2e3)
            sd_1e3 = float(sd_1e3)
            sd_5e2 = float(sd_5e2)
            sd_2e2 = float(sd_2e2)

            qpcr_values.append([qpcr_kit, gene_target, mean_2e3, sd_2e3, n_labs_2e3, mean_1e3, sd_1e3, n_labs_1e3, mean_5e2, sd_5e2, n_labs_5e2, mean_2e2, sd_2e2,n_labs_2e2])


    df = pd.DataFrame(qpcr_values, columns=["qpcr_kit", "gene_target", "mean_2e3", "sd_2e3", "n_labs_2e3", "mean_1e3", "sd_1e3", "n_labs_1e3", "mean_5e2", "sd_5e2", "n_labs_5e2", "mean_2e2", "sd_2e2", "n_labs_2e2"])

    # groupby gene_target. Then create the means for each mean and standard deviation
    print(df.groupby("gene_target")[["mean_2e3", "sd_2e3", "mean_1e3", "sd_1e3", "mean_5e2", "sd_5e2", "mean_2e2", "sd_2e2"]].mean())

