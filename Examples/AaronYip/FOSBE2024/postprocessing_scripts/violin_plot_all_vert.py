import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.rcParams.update({'font.size': 16})

# Import csv
path_to_csv = r'C:\Users\aaron\OneDrive - University of Waterloo\Research\Papers\Biophysics Calibration\simulation_data\abc_validation'
csv_file = 'all_summary_stats.csv'
df = pd.read_csv(os.path.join(path_to_csv, csv_file))

# Melt the DataFrame
DataRS = pd.melt(df, id_vars=["Exp", "gamma", "reg_param"])

# Plot
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(7, 7))
sns.violinplot(data=DataRS, x='variable', y='value', hue='Exp', inner='box', dodge=True, linewidth=1.5, palette='colorblind', scale='width', ax=ax)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.11), ncol=3)
sns.stripplot(data=DataRS, x='variable', y='value', hue='Exp', dodge=True, linewidth=0.5, alpha=0.5, jitter=True, palette='colorblind', size=4, legend=False)

# Separate summary stats with vertical lines
vertlines = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
for x in vertlines:
    plt.axvline(x = x, color = 'grey', alpha=0.5)

plt.ylim([-0.05,1.12])
ax.set_xticklabels(['Aspect Ratio', 'Convexity', 'Order', 'Density', 'Exp. Fit', 'Growth Rate'])
plt.yticks([0, 0.25, 0.50, 0.75, 1.00])
plt.xticks(rotation=45)
plt.ylabel('')
plt.xlabel('')
plt.tight_layout()
plt.show()
