import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams.update({'font.size': 11})

# Import csv
path_to_csv = r'C:\Users\aaron\OneDrive - University of Waterloo\Research\Papers\Biophysics Calibration\simulation_data\abc_validation'
csv_file = 'summary_stats_calibrated.csv'
df = pd.read_csv(os.path.join(path_to_csv, csv_file))

# Get data
data = df[['aspect_ratio', 'convexity', 'order_parameter', 'density_parameter', 'R2_exp_fit']]

# Create a figure and axes
fig, ax = plt.subplots(figsize=(4,4))

# Create a violin plot
quantiles = [0,0.25,0.5,0.75,1]
quantile_list = [quantiles for i in data.columns]
print(quantiles)
violin_parts = ax.violinplot(data, showmedians=False, quantiles=quantile_list)

# Add scatter points on top of the violin plot
for i, col in enumerate(data.columns):
    y = np.random.normal(i + 1, 0.04, size=len(data))
    ax.scatter(y, data[col], s=15, color='white', edgecolor='black', zorder=10, alpha=1.0)

# Customize the plot
plt.ylim(bottom=0, top=1.05)
ax.set_ylabel('')
ax.set_xticks(np.arange(1, len(data.columns) + 1))
ax.set_xticklabels(data.columns, rotation=90)
for pc in violin_parts['bodies']:
    pc.set_facecolor('#D43F3A')
    pc.set_edgecolor('black')
    pc.set_alpha(1)
    
cquantile_colors = violin_parts['cquantiles'].get_color()
colors = ['black' for i in quantiles]
violin_parts['cquantiles'].set_color(colors)

cbar_colors = violin_parts['cbars'].get_color()
colors = ['black' for i in data.columns]
violin_parts['cbars'].set_color(colors)

# Add a legend for scatter points
scatter_legend = ax.scatter([], [], color='white', edgecolor='black', label='Scatter Points')
#ax.legend(handles=[scatter_legend], loc='lower right')

# Label summary stats
ax.set_xticklabels(['Aspect Ratio', 'Convexity', 'Order', 'Density', 'Exp. Fit'])
plt.xticks(rotation=45)

# Show the plot
plt.tight_layout()
plt.show()
