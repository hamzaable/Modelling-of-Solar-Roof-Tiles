# -*- coding: utf-8 -*-
"""
Created on Thu May 12 12:19:06 2022

@author: mariu
"""

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

plt.close('all')


# Abbildung 5
d = {'electric energy gain (%)': [3.34, 3.33, 3.29, 3.22, 3.11, 0], 'AEP electric AC (kWh/a)': [2921, 2921, 2920, 2918, 2915, 2827]}
df = pd.DataFrame(data=d)

matplotlib.rc_file_defaults()
ax1 = sns.set_style(style="whitegrid", rc=None )


fig, ax1 = plt.subplots(figsize=(8,5))
sns.barplot(data = df, x=['$30 m^3/h$', '$25 m^3/h$', '$20 m^3/h$', '$15 m^3/h$', '$10 m^3/h$', '-'], y='AEP electric AC (kWh/a)', palette="Blues_d", alpha=0.5,  ax=ax1)
ax1.bar_label(ax1.containers[0])

ax2 = ax1.twinx()
ax2.grid(False)


sns.lineplot(data = df['electric energy gain (%)'], marker='o', sort = False, ax=ax2, label="Additional electric energy gain (AC) - cooling effect")


ax2.legend(loc="lower center")
ax2.legend(bbox_to_anchor=(0.83, -0.08))

ax1.set_ylabel('AEP electric AC (kWh/a)',fontdict= { 'fontsize': 12.5}) # 'fontweight':'bold'
ax2.set_ylabel('electric energy gain (%)',fontdict= { 'fontsize': 12.5}) # 'fontweight':'bold'
plt.subplots_adjust(bottom=0.18)



# Abbildung 6 
pr_cl =[87.96, 88.91, 85.11, 86.46, 83.22, 82.42, 82.98, 83.29, 83.73, 86.50, 85.31, 87.90]
pr = [86.92, 87.08, 82.93, 85.09, 81.27, 81.13, 81.23, 81.10, 81.94, 84.48, 83.81, 87.18]

d = {'PR with solar cooling': pr_cl, 'PR': pr}
df = pd.DataFrame(data=d)
month = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'Oktober', 'November', 'December']


matplotlib.rc_file_defaults()

ax1 = sns.set_style(style="whitegrid", rc=None )

fig, ax1 = plt.subplots(figsize=(8,3.8))

#sns.lineplot(data = df, x=month, y=df['PR with solar cooling'], marker='^', linestyle="dashdot")
#ax1.set_ylabel('Performance Ratio (%)',fontdict= { 'fontsize': 12.5}) # 'fontweight':'bold', 

sns.lineplot(data=df, markers=['^', '^'],  dashes=[(2, 2), (2, 2)],  linewidth=3.0, markersize=12) #linestyles="dashdot",scale=2,sizes=[200,200]
ax1.set(ylim=(76, 90))
ax1.set_xticklabels(month, rotation=45, fontdict= { 'fontsize': 13})
#ax1.set(xticks=range(1, 12), xticklabels=month)
#plt.xticks(rotation=45)
ax1.set_xlabel('2021', fontdict= { 'fontsize': 13})
ax1.set_ylabel('Performance Ratio (%)', fontdict= { 'fontsize': 13})
ax1.set_yticklabels(ax1.get_yticks(), size = 12.5)
ax1.xaxis.grid(False)
plt.legend(fontsize=12.5,  markerscale=1.5)
#plt.plot(markerscale=1.5)
ax1.set(xticks=range(0, 12), xticklabels=month)
plt.subplots_adjust(bottom=0.325)
