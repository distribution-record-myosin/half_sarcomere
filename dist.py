import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.clf()

df = pd.read_csv(f"dist.csv", delim_whitespace=True, header=None)
x = df.iloc[:, 1]
y = df.iloc[:, 2]
#plt.ylim(-19,23)
#plt.xlim(-7,14)
plt.scatter(x, y,marker=".",s=1)
plt.xlabel("$x_L$",fontsize=20)
plt.ylabel("$x_S$",fontsize=20)
plt.title(f"states of bound myosins",fontsize=20)
plt.grid()
plt.tick_params(axis='x', labelsize=20)  # x軸の目盛りラベル
plt.tick_params(axis='y', labelsize=20)  # y軸の目盛りラベル
plt.tight_layout() 
plt.savefig("dist.png") 