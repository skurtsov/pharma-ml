#EDA for activities csv
from pandas import read_csv
import numpy as np
import matplotlib.pyplot as plt

df = read_csv("data/egfr_activities.csv")
df.info()
print(df["ic50"].describe())
#MATH
log_ic50 = np.log10(df["ic50"])

log_ic50.hist(bins=50)

plt.xlabel("log10(IC50 nM)")
plt.ylabel("Number of molecules")
plt.title("EGFR IC50 Distribution")
plt.axvline(
    x=3,
    linestyle="--",
    label="Activity threshold (1000 nM)",
    color="red"
)

plt.legend()
plt.show()

log_ic50 = np.log10(df["ic50"])
