
from matplotlib import pyplot as plt
import numpy as np
import argparse

from sklearn import linear_model

parser = argparse.ArgumentParser(description="Arrhenous plot")

parser.add_argument("path", default="time.csv", help="Path to csv formatted input file")

args = parser.parse_args()

fname = args.path

temp = []
time = []

with open(fname) as file:
    for line in file:
        bits = line.split()
        time.append(float(bits[-1]))
        temp.append(float(bits[-2]))

temp = np.asarray(temp)
time = np.asarray(time)

plt.semilogy(1 / temp, time, "k+", label="Data")

x = 1 / temp
y = np.log(time)

# t = t0 * exp(dE / kT)
# ln(t) = ln(t0) + dE / kT

# ---------- OLS regression ----------

z, cov = np.polyfit(x, y, 1, cov = True)
stddev = np.sqrt(np.diag(cov))

kB = 8.617333262e-5
dE = abs(z[0]) * kB
err = stddev[0] * kB

print(f"   OLS: dE = {dE} pm {err}")

plt.semilogy(1 / temp, np.exp(np.polyval(z, x)), label="OLS")

# ---------- Robust ----------


xp = x.reshape(-1, 1)
yp = y.reshape(-1, 1)
trials = 100

best = 0
y_ransac = None

dEs = []

for _ in range(trials):

    ransac = linear_model.RANSACRegressor(max_trials=500)
    ransac.fit(xp, yp)
    
    dE = abs(ransac.estimator_.coef_[0][0]) * kB
    dEs.append(dE)

    score = ransac.score(xp, yp)

    if score > best:
        best = score
        y_ransac = ransac.predict(xp)

dE = np.mean(dEs)
err = np.std(dEs)

print(f"RANSAC: dE = {dE} pm {err}")

plt.semilogy(1 / temp, np.exp(y_ransac), label="RANSAC")



# ---------- IO ----------

plt.legend()
plt.savefig("arr.pdf")



# plt.show()



