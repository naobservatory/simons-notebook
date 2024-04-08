n_swabs = [1000, 500, 200, 100, 50]
studies = ["Lu et al. 2021", "Babiker et al. 2020", "Mostafa et al. 2020"]

violin_df = pd.DataFrame()

for i, (study, positive_ras) in enumerate(
    zip(studies, [df_op_lu_ras, df_np_babiker_ras, df_np_mostafa_ras])
):
    study_df = simulate_simple(positive_ras, n_swabs)
    study_df_log = np.log10(study_df)
    study_df_log_melted = study_df_log.reset_index().melt(
        id_vars=["index"],
        var_name="Number of Swabs",
        value_name="Log10 Relative Abundance",
    )
    study_df_log_melted["Study"] = study
    violin_df = pd.concat([violin_df, study_df_log_melted], ignore_index=True)


violin_df_100_swabs = violin_df[violin_df["Number of Swabs"] == 100]


fig, ax = plt.subplots(figsize=(8, 3))
sns_default_colors = sns.color_palette()
sns.violinplot(
    ax=ax,
    data=violin_df_100_swabs,
    x="Log10 Relative Abundance",
    y="Study",
    bw=0.5,
    palette="viridis",
    scale="width",
    inner=None,
    width=0.5,
)
ax.set_ylabel("")
ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_func))
for x in np.arange(-10, 0, 2):
    ax.axvline(x=x, color="black", linestyle="--", linewidth=0.5, alpha=0.5, zorder=-1)
