import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.DataFrame()

df['x'] = [250, 500, 1000, 1500]
df['RAxML-NG'] = [0.358970396742, 0.281401328493, 0.206761354298, 0.188365444752]
df['Treerecs'] = [0.258528783379, 0.187579630999, 0.129546138407, 0.119346474143]
df['Phyldog'] = [0.272231488512, 0.20714293281, 0.131212412268, 0.131669602346]
df['Notung'] = [0.355859146335, 0.218086125849, 0.168773155905, 0.157012245646]
df['JointSearch'] = [0.20784454721, 0.165717366227, 0.0997375796799, 0.0911058132789]

print(df)

# Plot the responses for different events and regions
ax = sns.lineplot(   data=df)
ax.set(xlabel="sites", ylabel="average RF distance")


plt.show()

