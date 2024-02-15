
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
msd = []
msd_c = []

with open(fname) as file:
    for line in file:
        bits = line.split()
        time.append(float(bits[-1]))
        temp.append(float(bits[-2]))
        msd_c.append(float(bits[-3]))
        msd.append(float(bits[-4]))

temp = np.asarray(temp)
time = np.asarray(time)
msd = np.asarray(msd)
msd_c = np.asarray(msd_c)

diff = msd_c / (6 * time)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.semilogy(1 / temp, time, "r+", label="Life")
ax1.set_xlabel("1 / T")
ax1.set_ylabel("Lifetime")

ax2.semilogy(1 / temp, diff, "b+", label="D_eff")
ax2.set_ylabel("D_eff")

ax1.set_title(fname)


# ---------- Regression

kB = 8.617333262e-5


x = 1 / temp
y_t = np.log(time)
y_d = np.log(diff)


def fit_ols(x, y):
    z, cov = np.polyfit(x, y, 1, cov = True)
    stddev = np.sqrt(np.diag(cov))
    dE = abs(z[0]) * kB
    err = stddev[0] * kB
    return dE, err, np.polyval(z, x)


dE, err, y_t_ols = fit_ols(x, y_t)
print(f"life  OLS: dE = {dE} pm {err}")
ax1.semilogy(1 / temp, np.exp(y_t_ols), "r-", label="Life - ols")

dE, err, y_d_ols = fit_ols(x, y_d)
print(f"diff  OLS: dE = {dE} pm {err}")
ax2.semilogy(1 / temp, np.exp(y_d_ols), "b-", label="D_eff - ols")


# ---------- Robust regression ----------


def fit_robust(x, y, trials = 100):

    xp = x.reshape(-1, 1)
    yp = y.reshape(-1, 1)

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

    return dE, err, y_ransac

dE, err, y_t_robust = fit_robust(x, y_t)
print(f"life RBST: dE = {dE} pm {err}")
ax1.semilogy(1 / temp, np.exp(y_t_robust), "k--", label="Robust")

dE, err, y_d_robust = fit_robust(x, y_d)
print(f"diff RBST: dE = {dE} pm {err}")
ax2.semilogy(1 / temp, np.exp(y_d_robust), "k--")

# ---------- IO ----------

fig.subplots_adjust(top=0.85)
fig.legend(ncol=5, loc="upper center")
fig.savefig("dis.pdf")

exit(0)

# ---------- IO ----------


# plt.show()



