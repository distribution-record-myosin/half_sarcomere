import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df1 = pd.read_csv(f"trans.csv", delim_whitespace=True, header=None)
x_arr1=df1.iloc[:,8].astype(float).tolist()
x_arr1 = np.array(x_arr1) 

plt.plot(x_arr1,linestyle='None',marker='o',ms=1)
plt.xlabel("# of 0.1 ms")
plt.grid()
plt.ylabel("Half-sarcomea shortening length [nm]")
plt.savefig("hssl.png") 
