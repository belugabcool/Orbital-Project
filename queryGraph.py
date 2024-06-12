import pandas as pd
import matplotlib.pyplot as plt

ind = 1

fig, axes = plt.subplots(nrows=2, ncols=2)

df = pd.read_csv(f"data\\dat{ind}.csv")
df.plot(x="t", y="V", ax=axes[0,0])
df.plot(x="X", y="V", ax=axes[0,1])
df.plot(x="t", y="theta", ax=axes[1,0])
df.plot(x="X", y="theta", ax=axes[1,1])


plt.show()