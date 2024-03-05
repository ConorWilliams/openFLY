
from matplotlib import pyplot as plt
import numpy as np
import argparse
from sys import float_info

from sklearn import linear_model

parser = argparse.ArgumentParser(description="Arrhenous plot")

parser.add_argument("path", default="time.csv", help="Path to csv formatted input file")
parser.add_argument("-o", "--output", default="plot.pdf", help="Output file")
parser.add_argument("-m", "--msd", action="store_true", help="Plot msd instead of lifetime")

args = parser.parse_args()

fname = args.path

temp = []
time = []
msd_h = []
msd_f = []

with open(fname) as file:
    for line in file:
        bits = line.split()

        if float(float(bits[-2])) < 470:
            continue

        time.append(float(bits[-1]))
        temp.append(float(bits[-2]))

        if (args.msd):
            msd_f.append(float(bits[-3]) * 1e-20)  # Angstroms^2
            msd_h.append(float(bits[-4]) * 1e-20)    # Angstroms^2

temp = np.asarray(temp)
inv_temp = 1 / temp
time = np.asarray(time)

if (args.msd):
    msd_h = np.asarray(msd_h)
    msd_f = np.asarray(msd_f)
    diff = msd_f / (6 * time)

fig, ax1 = plt.subplots()

if (args.msd):

    ax1.yaxis.label.set_color("red")
    ax1.tick_params(axis="y", labelcolor="red")

    ax2 = ax1.twinx()
    ax2.yaxis.label.set_color("blue")
    ax2.tick_params(axis="y", labelcolor="blue")

ax1.semilogy(inv_temp, time, "r+", label="Lifetime")
ax1.set_xlabel("Inverse-temperature/K^-1")
ax1.set_ylabel("Lifetime")

ax1.set_ylim(1e-11, 1e-5)

if (args.msd):
    ax2.semilogy(inv_temp, diff, "b+", label="Diffusivity")
    ax2.set_ylabel("D_eff")

# ---------- Second x - axis 

def fix_zero(some_array):
    corrected_array = some_array.copy()
    for i, entry in enumerate(some_array):
        # If element is zero, set to some small value
        if abs(entry) < float_info.epsilon:
            corrected_array[i] = float_info.epsilon
    
    return corrected_array

def invert(x):
    return 1 / fix_zero(x)

ax1.invert_xaxis()

axT = ax1.secondary_xaxis('top', functions=(invert, invert))
axT.set_xlabel('Temprature/K', fontsize='large')

# ---------- Regression

kB = 8.617333262e-5


x = inv_temp
y_t = np.log(time)

if (args.msd):
    y_d = np.log(diff)


# ---------- OLS ----------

def fit_ols(x, y):
    z, cov = np.polyfit(x, y, 1, cov = True)
    return np.polyval, z, np.sqrt(np.diag(cov))

def bands(f, z, dz, x, n_std):
    
    lims = [
        f((z[0] + n_std * dz[0], z[1] + n_std * dz[1]), x),
        f((z[0] + n_std * dz[0], z[1] - n_std * dz[1]), x),
        f((z[0] - n_std * dz[0], z[1] + n_std * dz[1]), x),
        f((z[0] - n_std * dz[0], z[1] - n_std * dz[1]), x)
    ]

    np.stack(lims, axis=0)

    y_lo = np.min(lims, axis=0)
    y_hi = np.max(lims, axis=0)

    return y_lo, y_hi

def plot_area(ax, x, f, z, dz, color, label):
    # Remove duplicates
    x = np.unique(np.sort(x)) 
    ax.semilogy(x, np.exp(f(z, x)), "k--")
    y_lo, y_hi = bands(f, z, dz, x, 1) 
    ax.fill_between(x, np.exp(y_lo), np.exp(y_hi), color=color, alpha=0.2)

f, z, dz = fit_ols(x, y_t)
print(f"life  OLS: dE = {np.abs(z[0]) * kB:.4f} pm {dz[0] * kB:.4f}")
plot_area(ax1, x, f, z, dz, "r", "Lifetime")


if (args.msd):
    f, z, dz = fit_ols(x, y_d)
    i300 = 1 / 300.
    D300 = np.exp(f(z, i300))
    D300_range = np.exp(bands(f, z, dz, i300, 1))
    D300_err = (D300_range[1] - D300_range[0]) / 2

    print(f"diff  OLS: dE = {np.abs(z[0]) * kB:.4f} pm {dz[0] * kB:.4f}, D_eff(300)={D300:.3e} pm {D300_err:.3e}")    
    plot_area(ax2, x, f, z, dz, "b", "Diffusivity")


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
        
        dEs.append(abs(ransac.estimator_.coef_[0][0]) * kB)

        score = ransac.score(xp, yp)

        if score > best:
            best = score
            y_ransac = ransac.predict(xp)

    dE = np.mean(dEs)
    err = np.std(dEs)

    return dE, err, y_ransac

dE, err, y_t_robust = fit_robust(x, y_t)
print(f"life RBST: dE = {dE:.4f} pm {err:.4f}")
# ax1.semilogy(x, np.exp(y_t_robust), "k--", label="Robust")

if (args.msd):
    dE, err, y_d_robust = fit_robust(x, y_d)
    print(f"diff RBST: dE = {dE:.4f} pm {err:.4f}")
    # ax2.semilogy(x, np.exp(y_d_robust), "k--")

# ---------- IO ----------

fig.subplots_adjust(right=0.875)
fig.savefig(args.output, format="pdf")

exit(0)

# ---------- IO ----------


# plt.show()



