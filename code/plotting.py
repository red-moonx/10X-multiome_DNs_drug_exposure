import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_gene_in_time(
    adata, 
    gene, 
    time_key='timepoint', 
    batch_key='replicate',
    save_path=None,
    save=False,
):
    # Define colors for each replicate number
    replicate_colors = {
        'R1': "#8B2635",
        'R2': "#3E9EB6",
        'R3': "#729C44",
        'R4': "#BC8034",
    }

    # initiate plot
    fig, axs = plt.subplots(1, 1, figsize=(8, 3))

    timepoints = adata.obs[time_key].cat.categories
    std_per_timepoint = []
    mean_per_timepoint = []

    # get mean and std per timepoint
    for time in timepoints:
        std = np.std([
            adata[(adata.obs[time_key] == time) &
                  (adata.obs[batch_key] == i)][:, gene].X.mean()
            for i in adata[adata.obs[time_key] == time].obs[batch_key].cat.categories])
        std_per_timepoint.append(std)
        mean = adata[(adata.obs[time_key] == time)][:, gene].X.mean()
        mean_per_timepoint.append(mean)

    # plot mean and std (as errorbar) in gray with alpha
    axs.errorbar(timepoints, mean_per_timepoint, yerr=std_per_timepoint, capsize=6, color="gray", alpha=0.6)
    axs.plot(timepoints, mean_per_timepoint, color="gray", alpha=0.6)

    # plot each replicate as dot with specific color
    for time in timepoints:
        replicates = adata[adata.obs[time_key] == time].obs[batch_key].cat.categories
        for replicate in replicates:
            replicate_number = replicate.split(' ')[-1]  # Extract replicate number from the string
            replicate_color = replicate_colors.get(replicate_number, "red")  # Default to red if replicate_number is not found
            axs.scatter([time], [adata[(adata.obs[time_key] == time) & (adata.obs[batch_key] == replicate)][:, gene].X.mean()],
                        c=replicate_color)

    axs.set_title(gene)
    axs.grid(False)
    axs.yaxis.tick_right()

    # remove spines
    axs.spines['left'].set_visible(False)
    axs.spines['top'].set_visible(False)
    
    # Remove axis labels
    fig.texts.clear()
    
    # Remove X-axis tick labels but keep the ticks
    axs.set_xticklabels([])

    # save figure
    if save:
        fig.savefig(save_path, bbox_inches='tight', format='pdf', dpi=300)

    plt.show()