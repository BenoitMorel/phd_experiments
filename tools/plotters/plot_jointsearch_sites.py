import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.DataFrame()

df['x'] = [250, 500, 1000]
df['RAxML-NG'] = [0.139209604236, 0.139209604236, 0.139209604236]
df['Treerecs'] = [0.0386519664884, 0.0386519664884, 0.0386519664884]
df['Phyldog'] = [0.0244567815894, 0.0244567815894, 0.0244567815894]
df['JointSearch'] = [0.017849850674, 0.0178498506740, 0.017849850674]

print(df)

# Plot the responses for different events and regions
ax = sns.lineplot(   data=df)
ax.set(xlabel="sites", ylabel="average RF distance")


plt.show()

