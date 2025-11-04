import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.rcParams.update({'font.size': 14})

# Import csv
path_to_csv = r'C:\Users\aaron\OneDrive - University of Waterloo\Research\Papers\Biophysics Calibration\simulation_data\abc_validation'
csv_file = 'all_summary_stats.csv'
df = pd.read_csv(os.path.join(path_to_csv, csv_file))

# Melt the DataFrame
DataRS = pd.melt(df, id_vars=["Exp", "gamma", "reg_param"])

# Plot
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(8, 8))
sns.violinplot(data=DataRS, x='value', y='variable', hue='Exp', inner='box', dodge=True, linewidth=1.5, palette='colorblind', scale='width', ax=ax)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3)

sns.stripplot(data=DataRS, x='value', y='variable', hue='Exp', dodge=True, linewidth=0.5, alpha=0.5, jitter=True, palette='colorblind', size=4, legend=False)

plt.xlim([-0.05,1.12])
#ax.grid(False)
ax.set_yticklabels(['Aspect Ratio', 'Convexity', 'Order', 'Density', 'Exp. Fit', 'Growth Rate'])
plt.xticks([0, 0.25, 0.50, 0.75, 1.00])
plt.yticks(rotation=45)
plt.xlabel('')
plt.show()
