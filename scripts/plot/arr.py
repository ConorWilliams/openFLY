
from matplotlib import pyplot as plt
import numpy as np
import argparse

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

# t = t0 * exp(dE / kT)
# ln(t) = ln(t0) + dE / kT

x = 1 / temp
y = np.log(time)
z, cov = np.polyfit(x, y, 1, cov = True)

stddev = np.sqrt(np.diag(cov))

kB = 8.617333262e-5

dE = abs(z[0]) * kB
err = stddev[0] * kB



print(f"dE = {dE} pm {err}")
# print(f"dE = {dE-err}-{dE+err}")

plt.semilogy(1 / temp, time, "k+")
plt.semilogy(1 / temp, np.exp(np.polyval(z, x)))


plt.savefig("arr.pdf")



# plt.show()



