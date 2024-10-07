"""
EDAM : Element Definition AutoMator
Yaël MOUSSOUNI
University of Strasbourg, Department of Physics and Engineering
"""
# Imports
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
# Variables and parameters
var_X = "PHI"
prefix = "EDAM_simple_{}_".format(var_X) # Prefix of the cfg_files
super_prefix = "EDAM_ncal_O4_{}_".format(var_X)
directory = "cfg_files/EDAM/" # working directory for cfg files
out_dir = "./"

sep = ";" # CSV columns
endl = "\n" # CSV lines

superposition = True
analytical = True
N = 1000 # nb of points for analytical model

r, phi, z = 1, 0, 0 # default parameters
# m, DEG, m
m = 2.288 # kg
M = 42.3648 # kg
R_prime = 66.217e-3 # m - rotor masse center distance from the axis

to_png = True # export to png ?
to_screen = True # show the result ?
# File reading
file = open(out_dir + prefix+"output.csv", "r")
line1 = file.readline()
lines = file.readlines()
file.close()

if superposition :
   super_file = open(out_dir + super_prefix+"output.csv", "r")
   super_line1 = super_file.readline()
   super_lines = super_file.readlines()
   super_file.close()
# Labels
labels = line1.replace(endl, "").split(sep)
xlabel = labels[0]
y1label = labels[1]
y2label = labels[2]
y3label = labels[3]
plot1label = "simple model"
plot2label = "more realistic model"
plot3label = "analytical model"
# Plot parameters
plt.style.use("default")
fig, axs = plt.subplots(3)
fig.set_size_inches(8.0, 7.0)
ax1 = axs[0]
ax2 = axs[1]
ax3 = axs[2]

plt.suptitle("EDAM, a FROMAGE extension")
ax1.set_title("Force offset with rotor position")
ax1.set_xlabel(xlabel)
ax1.set_ylabel(y1label)
ax1.get_yaxis().get_major_formatter().set_useOffset(False)
ax2.set_title("Force at $2f$ with rotor position")
ax2.set_xlabel(xlabel)
ax2.set_ylabel(y2label)
ax2.get_yaxis().get_major_formatter().set_useOffset(False)
ax3.set_title("Force amplitude with rotor position")
ax3.set_xlabel(xlabel)
ax3.set_ylabel(y3label)
ax3.get_yaxis().get_major_formatter().set_useOffset(False)

# ax1.set_yscale('symlog')
# ax2.set_yscale('symlog')
# ax3.set_yscale('symlog')
# Data extraction
X = []
F0 = []
F2 = []
A = []
for line in lines:
   x, f0, f2, a = line.replace(endl, "").split(sep)
   X.append(float(x))
   F0.append(float(f0))
   F2.append(float(f2))
   A.append(float(a))
X = np.array(X)
F0 = np.array(F0)
F2 = np.array(F2)
A = np.array(A)

if superposition:
   super_X = []
   super_F0 = []
   super_F2 = []
   super_A = []
   for super_line in super_lines:
      super_x, super_f0, super_f2, super_a = super_line.replace(endl, "").split(sep)
      super_X.append(float(super_x))
      super_F0.append(float(super_f0))
      super_F2.append(float(super_f2))
      super_A.append(float(super_a))
   super_X = np.array(super_X)
   super_F0 = np.array(super_F0)
   super_F2 = np.array(super_F2)
   super_A = np.array(super_A)
# analytical solution
def d0_definition(r0, z0):
   return np.sqrt(r0**2+z0**2)

def F_min(m, M, d_0, R_prime, phi, xi):
   return 2*const.G*m*M * 1 / (d_0**2 + R_prime**2) * np.cos(phi/180*np.pi) * np.cos(xi)

def F_max(m, M, d_0, R_prime, phi, xi):
   return  2*const.G*m*M * (d_0**2 + R_prime**2) / (d_0**2 - R_prime**2)**2 * np.cos(phi/180*np.pi)

def F0_analytique(m, M, d_0, R_prime, phi, xi):
   return (F_max(m, M, d_0, R_prime, phi, xi) + F_min(m, M, d_0, R_prime, phi, xi))/2

def F2_analytique(m, M, d_0, R_prime, phi, xi):
   return (F_max(m, M, d_0, R_prime, phi, xi) - F_min(m, M, d_0, R_prime, phi, xi))/2

# data plot
ax1.plot(X, F0, ".", color='#003763', markersize=4, label=plot1label)
ax2.plot(X, F2, ".", color='#003763', markersize=4, label=plot1label)
ax3.plot(X, A, ".", color='#003763', markersize=4, label=plot1label)
if superposition:
   ax1.plot(super_X, super_F0, "o", color='#AA10AA', markersize=5, label=plot2label)
   ax2.plot(super_X, super_F2, "o", color='#AA10AA', markersize=5, label=plot2label)
   ax3.plot(super_X, super_A, "o", color='#AA10AA', markersize=5, label=plot2label)
if analytical:
   x_min = np.min(X)
   x_max = np.max(X)
   a_X = np.linspace(x_min, x_max, N)
   if var_X == "R":
      R = a_X
      D0 = d0_definition(R, z)
      XI_MIN = np.pi/4 - np.arctan((D0-R_prime)/(D0+R_prime))
      a_F0 = F0_analytique(m, M, D0, R_prime, phi, XI_MIN)
      a_F2 = F2_analytique(m, M, D0, R_prime, phi, XI_MIN)
   if var_X == "PHI":
      PHI = a_X
      d0 = d0_definition(r, z)
      xi_min = np.pi/4 - np.arctan((d0-R_prime)/(d0+R_prime))
      a_F0 = F0_analytique(m, M, d0, R_prime, PHI, xi_min)
      a_F2 = F2_analytique(m, M, d0, R_prime, PHI, xi_min)
   if var_X == "Z":
      Z = a_X
      D0 = d0_definition(r, Z)
      XI_MIN = np.pi/4 - np.arctan((D0-R_prime)/(D0+R_prime))
      a_F0 = F0_analytique(m, M, D0, R_prime, phi, XI_MIN)
      a_F2 = F2_analytique(m, M, D0, R_prime, phi, XI_MIN)
   a_A = a_F2
   ax1.plot(a_X, a_F0, "-", color='#10AA10', label=plot3label)
   ax2.plot(a_X, a_F2, "-", color='#10AA10', label=plot3label)
   ax3.plot(a_X, a_A, "-", color='#10AA10', label=plot3label)

ax1.legend()
ax2.legend()
ax3.legend()
ax1.grid()
ax2.grid()
ax3.grid()
plt.tight_layout()
if to_png : plt.savefig("edam")
if to_screen : plt.show()