import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns

# https://python-graph-gallery.com/122-multiple-lines-chart/

df = pd.DataFrame()
f, ax = plt.subplots(1)
ax.set_ylim(ymin=0)
methods = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "JointSearch"]

df['bl'] = [0.5, 1.0, 2.0, 4.0]
df['RAxML-NG'] = [0.139552240776, 0.341051205447, 0.572642826728, 0.680660436994]
df['Treerecs'] = [0.0731537739958, 0.15988898555, 0.280588063666, 0.404015266016]
df['Phyldog'] = [0.0876450176615, 0.222595175173, 0.432981610856, 0.540814984988]
df['Notung'] = [0.0944926568377, 0.221493198255, 0.422827046051, 0.514044710917]
df['JointSearch'] = [0.0548871148658, 0.139352516076, 0.251452327924, 0.357417753002]

for method in methods:
  plt.plot( 'bl', method, data=df, marker='x', linewidth=2)

plt.xlabel('Average BL factor')
plt.ylabel('RF distance')
plt.title("24 species, sites = 500, D=0.5 L=0.25")
plt.legend()

plt.show()

