"""
EDAM : Element Definition AutoMator
Yaël MOUSSOUNI
University of Strasbourg, Department of Physics and Engineering
"""
# Imports
import os
import numpy as np
import subprocess

# Variables and parameters
type_model = "simple" # available : simple, ncal_O4
var_X = "Z" # available : R, PHI, Z
prefix = "EDAM_{}_{}_".format(type_model, var_X) # Prefix of the cfg_files
directory = "cfg_files/EDAM/" # working directory for cfg files
out_dir = "./"

output_csv = True # output a CSV ?
sep = ";" # CSV columns
endl = "\n" # CSV lines
r, phi, z = 1, 0, 0
X_min = 0 # minimum value of X
X_max = 4 # maximum value of X
N = 100 # number of values of X

## FROMAGE
STEP = "10 36" # [angle] [number of steps] 
ARM_LENGTH = "3000"
SIGNAL = "2"

# File type
def simple_model(r, phi, z):
   name = "simple_model"
   ROTOR_CYLINDRICAL = str(r) + " " + str(phi) + " " + str(z)
   data = """# EDAM config file
### MIRROR DEFINITION
GRID_SIZE 1 1 1
CUBOID 42.3648 1 1 1 0 0 0
### ROTOR DEFINITION
ROTOR_CYLINDRICAL {} 0
GRID_SIZE 1 1 1
CUBOID 2.288 1 1 1 0 0.05702003335661673 0
CUBOID 2.288 1 1 1 0 -0.05702003335661673 0

### GENERAL PARAMETERS
STEP {}
ARM_LENGTH {}
SIGNAL {}
""".format(ROTOR_CYLINDRICAL, STEP, ARM_LENGTH, SIGNAL)
   return name, data

def ncal_O4_model(r, phi, z):
   name = "ncal_O4_model"
   ROTOR_CYLINDRICAL = str(r) + " " + str(phi) + " " + str(z)
   data = """# EDAM config file
### MIRROR DEFINITION
GRID_SIZE 12 30 8
CYLINDER 2202. 0 0.175 0.2 360 0 0 0

### ROTOR DEFINITION
ROTOR_CYLINDRICAL {} 0
GRID_SIZE 8 65 40
CYLINDER 2808.1 0.02900 0.104000 0.104400 90 0 0 0
CYLINDER 2808.1 0.02900 0.104000 0.104400 90 0 0 180

### GENERAL PARAMETERS
STEP {}
ARM_LENGTH {}
SIGNAL {}
""".format(ROTOR_CYLINDRICAL, STEP, ARM_LENGTH, SIGNAL)
   return name, data

# File generation
def file_gen(name, data):
   input_file = open(directory + prefix + name + ".cfg", "w")
   input_file.write(data)
   input_file.close()
   return None

# Type of variation
def r_var(model, phi, z, R):
   N = len(R)
   F0 = np.zeros(N)
   F2 = np.zeros(N)
   A = np.zeros(N)
   for i in range(N):
      r = R[i]
      name, data = model(r, phi, z)
      name += "_" + str(i)
      file_gen(name, data)
      fro_out = FROMAGE(directory, name)
      f0, f2, a = force_extract(fro_out)
      F0[i] = f0
      F2[i] = f2
      A[i] = f2
   return F0, F2, A

def phi_var(model, r, z, PHI):
   N = len(PHI)
   F0 = np.zeros(N)
   F2 = np.zeros(N)
   A = np.zeros(N)
   for i in range(N):
      phi = PHI[i]
      name, data = model(r, phi, z)
      name += "_" + str(i)
      file_gen(name, data)
      fro_out = FROMAGE(directory, name)
      f0, f2, a = force_extract(fro_out)
      F0[i] = f0
      F2[i] = f2
      A[i] = a
   return F0, F2, A

def z_var(model, r, phi, Z):
   N = len(Z)
   F0 = np.zeros(N)
   F2 = np.zeros(N)
   A = np.zeros(N)
   for i in range(N):
      z = Z[i]
      name, data = model(r, phi, z)
      name += "_" + str(i)
      file_gen(name, data)
      fro_out = FROMAGE(directory, name)
      f0, f2, a = force_extract(fro_out)
      F0[i] = f0
      F2[i] = f2
      A[i] = a
   return F0, F2, A
# FROMAGE execution
def FROMAGE(directory, name):
   fro_out = subprocess.run(["./fromage", "./{}".format(directory+prefix+name+".cfg")], capture_output=True).stdout.decode("utf-8") 
   return fro_out

# Data finder
def force_extract(fro_out):
   trig_amplitude = "Amplitude of the longitudinal force is "
   trig_offset = "Force offset is "
   trig_2f = "Force at 2f is "
   trig_unit = "N"
   i = fro_out.find(trig_amplitude) + len(trig_amplitude)
   j = fro_out[i:].find(trig_unit) + i
   k = fro_out.find(trig_offset) + len(trig_offset)
   l = fro_out[k:].find(trig_unit) + k
   m = fro_out.find(trig_2f) + len(trig_2f)
   n = fro_out[m:].find(trig_unit) + m
   A = float(fro_out[i:j-1])  # amplitude
   F0 = float(fro_out[k:l-1]) # offset
   F2 = float(fro_out[m:n-1]) # fourier
   return F0, F2, A

# CSV output
def csv_output(X, F0, F2, A, type_X="$X$"):
   output_file = open(out_dir + prefix+"output.csv", "w")
   output_file.write(type_X+sep+"Force offset $F_0$/N"+sep+"Force at $2f$ $F_2$/N"+sep+"Force amplitude $A$/N"+endl)
   for i in range(N):
      output_file.write(str(X[i])+sep+str(F0[i])+sep+str(F2[i])+sep+str(A[i])+endl)
   output_file.close()
   return None

X = np.linspace(X_min, X_max, N)
if type_model == "simple" :
   model = simple_model
elif type_model == "ncal_O4":
   model = ncal_O4_model

if var_X == "R":
   type_X = "$r$/m"
   F0, F2, A = r_var(model, phi, z, X)
elif var_X == "PHI": 
   type_X = "$\\phi$/deg"
   F0, F2, A = phi_var(model, r, z, X)
elif var_X == "Z": 
   type_X = "$z$/m"
   F0, F2, A = z_var(model, r, phi, X)

if output_csv : csv_output(X, F0, F2, A, type_X)