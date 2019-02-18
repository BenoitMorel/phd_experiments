import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# https://python-graph-gallery.com/122-multiple-lines-chart/

df = pd.DataFrame()
f, ax = plt.subplots(1)
ax.set_ylim(ymin=0)
methods = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "JointSearch"]

df['sites'] = [250, 500, 1000, 1500]
df['RAxML-NG'] = [0.358970396742, 0.281401328493, 0.206761354298, 0.188365444752]
df['Treerecs'] = [0.258528783379, 0.187579630999, 0.129546138407, 0.119346474143]
df['Phyldog'] = [0.272231488512, 0.20714293281, 0.131212412268, 0.131669602346]
df['Notung'] = [0.355859146335, 0.218086125849, 0.168773155905, 0.157012245646]
df['JointSearch'] = [0.20784454721, 0.165717366227, 0.0997375796799, 0.0911058132789]

for method in methods:
  plt.plot( 'sites', method, data=df, marker='x', linewidth=2)

plt.xlabel('Sites')
plt.ylabel('RF distance')
plt.title("24 species, bl factor = 1, D=1.0 L=0.5")
plt.legend()

plt.show()

