import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("1dlap.csv")
x    = data["x"].tolist()
t    = data["t"].tolist()

plt.plot(x, t)
plt.grid(True)
plt.show()