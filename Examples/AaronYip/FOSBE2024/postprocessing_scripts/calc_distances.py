import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import kl_divergence

def distance_calculation(sim_stats, exp_stats):
    """
    Calculate Euclidean distance between simulations and experimental observations
    
    @param  sim_stats   Dict of simulated summary statistics
    @param  exp_stats   Dict of experimental summary statistics
    @return distance    Euclidean distance between simuluations and experiment
    """
    #try:
    distance_dict = {}
    for x in sim_stats.keys():
        distance_dict[x] = np.absolute(kl_divergence.get_kl_divergence(exp_stats[x], sim_stats[x]))
    print("Distance calc succeeded")
    #except:
    #    print("Distance calculation failed")
    #    distance = 100 # Reject if there is an error

    return distance_dict

matplotlib.rcParams.update({'font.size': 14})

# Import csv
path_to_csv = r'C:\Users\aaron\OneDrive - University of Waterloo\Research\Papers\Biophysics Calibration\simulation_data\abc_validation'
csv_file = 'all_summary_stats.csv'
df = pd.read_csv(os.path.join(path_to_csv, csv_file))

summary_stats = ['aspect_ratio', 'convexity', 'order_parameter', 'density_parameter', 'R2_exp_fit', 'growth_rate']
exp_csv_file =  os.path.join(path_to_csv, 'summary_stats_experiment.csv')
calibrated_csv_file = os.path.join(path_to_csv, 'summary_stats_calibrated.csv')
prior_csv_file = os.path.join(path_to_csv, 'summary_stats_uniform.csv')

exp_data = pd.read_csv(exp_csv_file)
calibrated_data = pd.read_csv(calibrated_csv_file)
prior_data = pd.read_csv(prior_csv_file)

exp_summary_stats = {}
calibrated_summary_stats = {} 
prior_summary_stats = {}

for stat in summary_stats:
    exp_summary_stats[stat] = np.asarray(exp_data[stat])
    calibrated_summary_stats[stat] = np.asarray(calibrated_data[stat])
    prior_summary_stats[stat] = np.asarray(prior_data[stat])
    
exp_v_calibrated_dist = distance_calculation(calibrated_summary_stats, exp_summary_stats)
exp_v_prior_dist = distance_calculation(prior_summary_stats, exp_summary_stats)
diffs = {'Calibrated': list(exp_v_calibrated_dist.values()), 'Prior': list(exp_v_prior_dist.values())}

x = np.arange(len(summary_stats))  # the label locations
width = 0.25  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(figsize=(6,4))
colors = ['#C38820', '#158B6A']
i = 0
for attribute, measurement in diffs.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute, color=colors[i], edgecolor='black')
    multiplier += 1
    i += 1

ax.set_ylabel('|KL Divergence|')
ax.set_xticks(x + width, summary_stats)
ax.set_xticklabels(['Aspect Ratio', 'Convexity', 'Order', 'Density', 'Exp. Fit', 'Growth Rate'])
ax.legend(loc='upper center', ncols=2)
plt.xticks(rotation=45)
plt.ylim([0,2.5])
plt.tight_layout()
plt.show()

'''
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
'''