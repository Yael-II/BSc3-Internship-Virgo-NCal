"""
EDAM : Element Definition AutoMator
Yaël MOUSSOUNI
University of Strasbourg, Department of Physics and Engineering
"""
# Imports
import os
import numpy as np
import subprocess
import pandas as pd

dir_path = "./"
csv_path = "../Fichiers/"
name_rotors = "rotors.csv"
name_uncert_simple = "uncert_simple_empty.csv"
name_uncert_ncal_O4 = "uncert_ncal_O4_empty.csv"

type_model = "simple" # available : simple, ncal_O4
prefix = "EDAM_uncert_" # Prefix of the cfg_files
directory = "cfg_files/EDAM/" # working directory for cfg files

## FROMAGE
STEP = "22.5 16" # [angle] [number of steps] 
ARM_LENGTH = "3000"
SIGNAL = "2"

index_uncert = 1 # 0 = F0, 1 = F2, 2 = A

R_inn = 29e-3 # m
rho_mirror = 2202 # kg m-3
R_mirror = (350.320e-3)/2 # m
b_mirror = 199.250e-3 # m
dR_mirror = (0.01e-3)/2 # m
db_mirror = 0.01e-3 # m

# CSV files --> DataFrame

os.chdir(dir_path)
rotors_df = pd.read_csv(csv_path+name_rotors, delimiter=",", index_col="NCal").head(6)
uncert_df_simple = pd.read_csv(csv_path+name_uncert_simple, delimiter=",", index_col="NCal").head(6)
uncert_df_ncal_O4 = pd.read_csv(csv_path+name_uncert_ncal_O4, delimiter=",", index_col="NCal").head(6)
uncert_df_analyt =  pd.read_csv(csv_path+name_uncert_simple, delimiter=",", index_col="NCal").head(6)

# Computation of m and R_prime and their associated uncertainties

rotors_df["m"] = rotors_df["rho"]*np.pi*rotors_df["b_rot"]*(rotors_df["R_rot"]**2-R_inn**2)/4
rotors_df["dm"] = rotors_df["m"]*(rotors_df["drho"]/rotors_df["rho"] + rotors_df["db_rot"]/rotors_df["b_rot"]) + rotors_df["rho"] * rotors_df["b_rot"] * np.pi/2*rotors_df["R_rot"]*rotors_df["dR_rot"]

rotors_df["R_prime"] = 4*np.sqrt(2)/(3*np.pi)*(rotors_df["R_rot"]**3-R_inn**3)/(rotors_df["R_rot"]**2-R_inn**2)
rotors_df["dR_prime"] = rotors_df["dR_rot"]

# NCals cfg models

def simple_model(G, m, M, r, phi, z, R_prime):
   name = "simple_model"
   ROTOR_CYLINDRICAL = str(r) + " " + str(phi) + " " + str(z)
   data = """# EDAM config file
### MIRROR DEFINITION
GRID_SIZE 1 1 1
CUBOID {} 1 1 1 0 0 0
### ROTOR DEFINITION
ROTOR_CYLINDRICAL {} 0
GRID_SIZE 1 1 1
CUBOID {} 1 1 1 0 {} 0
CUBOID {} 1 1 1 0 -{} 0

### GENERAL PARAMETERS
G_GRAVITY {}
STEP {}
ARM_LENGTH {}
SIGNAL {}
""".format(M, ROTOR_CYLINDRICAL, m, R_prime, m, R_prime, G, STEP, ARM_LENGTH, SIGNAL)
   return name, data


def ncal_O4_model(G, rho_mirror, R_mirror, b_mirror, rho, b_rot, R_inn, R_rot, r, phi, z):
   name = "ncal_O4_model"
   ROTOR_CYLINDRICAL = str(r) + " " + str(phi) + " " + str(z)
   MIRROR = "{} 0 {} {} 360 0 0 0".format(rho_mirror, R_mirror, b_mirror)
   SECTOR_1 = "{} {} {} {} 90 0 0 0".format(rho, R_inn, R_rot, b_rot)
   SECTOR_2 = "{} {} {} {} 90 0 0 180".format(rho, R_inn, R_rot, b_rot)
   data = """# EDAM config file
### MIRROR DEFINITION
GRID_SIZE 12 30 8
CYLINDER {}

### ROTOR DEFINITION
ROTOR_CYLINDRICAL {} 0
GRID_SIZE 8 65 40
CYLINDER {}
CYLINDER {}

### GENERAL PARAMETERS
G_GRAVITY {}
STEP {}
ARM_LENGTH {}
SIGNAL {}
""".format(MIRROR, ROTOR_CYLINDRICAL, SECTOR_1, SECTOR_2, G, STEP, ARM_LENGTH, SIGNAL)
   return name, data

# File generation

def file_gen(name, data):
   input_file = open(directory + prefix + name + ".cfg", "w")
   input_file.write(data)
   input_file.close()
   return None

# FROMAGE

def FROMAGE(directory, name):
   fro_out = subprocess.run(["./fromage", "./{}".format(directory+prefix+name+".cfg")], capture_output=True).stdout.decode("utf-8") 
   return fro_out

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

def uncert_simple_file_gen(rotors_df):
   G = rotors_df["G"].to_numpy()
   m = rotors_df["m"].to_numpy()
   M = rotors_df["M"].to_numpy()
   r = rotors_df["r"].to_numpy()
   phi = rotors_df["phi"].to_numpy()
   z = rotors_df["z"].to_numpy()
   R_prime = rotors_df["R_prime"].to_numpy()
   dG = rotors_df["dG"].to_numpy()
   dm = rotors_df["dm"].to_numpy()
   dM = rotors_df["dM"].to_numpy()
   dr = rotors_df["dr"].to_numpy()
   dphi = rotors_df["dphi"].to_numpy()
   dz = rotors_df["dz"].to_numpy()
   dR_prime = rotors_df["dR_prime"].to_numpy()
   # ref
   for i in range(len(rotors_df)):
      name, data = simple_model(G[i], m[i], M[i], r[i], phi[i], z[i], R_prime[i])
      name = "ref_" + name + str(i)
      file_gen(name, data)
      F_ref = force_extract(FROMAGE(directory, name))

      # G
      name, data = simple_model(G[i]-dG[i], m[i], M[i], r[i], phi[i], z[i], R_prime[i])
      name = "G_min_" + name + str(i)
      file_gen(name, data)
      F_G_min = force_extract(FROMAGE(directory, name))


      name, data = simple_model(G[i]+dG[i], m[i], M[i], r[i], phi[i], z[i], R_prime[i])
      name = "G_max_" + name + str(i)
      file_gen(name, data)
      F_G_max = force_extract(FROMAGE(directory, name))

      # m
      name, data = simple_model(G[i], m[i]-dm[i], M[i], r[i], phi[i], z[i], R_prime[i])
      name = "m_min_" + name + str(i)
      file_gen(name, data)
      F_m_min = force_extract(FROMAGE(directory, name))

      name, data = simple_model(G[i], m[i]+dm[i], M[i], r[i], phi[i], z[i], R_prime[i])
      name = "m_max_" + name + str(i)
      file_gen(name, data)
      F_m_max = force_extract(FROMAGE(directory, name))

      # M
      name, data = simple_model(G[i], m[i], M[i]-dM[i], r[i], phi[i], z[i], R_prime[i])
      name = "M_min_" + name + str(i)
      file_gen(name, data)
      F_M_min = force_extract(FROMAGE(directory, name))

      name, data = simple_model(G[i], m[i], M[i]+dM[i], r[i], phi[i], z[i], R_prime[i])
      name = "M_max_" + name + str(i)
      file_gen(name, data)
      F_M_max = force_extract(FROMAGE(directory, name))

      # r
      name, data = simple_model(G[i], m[i], M[i], r[i]-dr[i], phi[i], z[i], R_prime[i])
      name = "r_min_" + name + str(i)
      file_gen(name, data)
      F_r_min = force_extract(FROMAGE(directory, name))

      name, data = simple_model(G[i], m[i], M[i], r[i]+dr[i], phi[i], z[i], R_prime[i])
      name = "r_max_" + name + str(i)
      file_gen(name, data)
      F_r_max = force_extract(FROMAGE(directory, name))

      # phi
      name, data = simple_model(G[i], m[i], M[i], r[i], phi[i]-dphi[i], z[i], R_prime[i])
      name = "phi_min_" + name + str(i)
      file_gen(name, data)
      F_phi_min = force_extract(FROMAGE(directory, name))

      name, data = simple_model(G[i], m[i], M[i], r[i], phi[i]+dphi[i], z[i], R_prime[i])
      name = "phi_max_" + name + str(i)
      file_gen(name, data)
      F_phi_max = force_extract(FROMAGE(directory, name))

      # z
      name, data = simple_model(G[i], m[i], M[i], r[i], phi[i], z[i]-dz[i], R_prime[i])
      name = "z_min_" + name + str(i)
      file_gen(name, data)
      F_z_min = force_extract(FROMAGE(directory, name))

      name, data = simple_model(G[i], m[i], M[i], r[i], phi[i], z[i]+dz[i], R_prime[i])
      name = "z_max_" + name + str(i)
      file_gen(name, data)
      F_z_max = force_extract(FROMAGE(directory, name))

      # R_prime
      name, data = simple_model(G[i], m[i], M[i], r[i], phi[i], z[i], R_prime[i]-dR_prime[i])
      name = "R_prime_min_" + name + str(i)
      file_gen(name, data)
      F_R_prime_min = force_extract(FROMAGE(directory, name))
      
      name, data = simple_model(G[i], m[i], M[i], r[i], phi[i], z[i], R_prime[i]+dR_prime[i])
      name = "R_prime_max_" + name + str(i)
      file_gen(name, data)
      F_R_prime_max = force_extract(FROMAGE(directory, name))

      uncert_df_simple["F0"][uncert_df_simple.index[i]] = F_ref[0]
      uncert_df_simple["F2"][uncert_df_simple.index[i]] = F_ref[1]
      uncert_df_simple["A"][uncert_df_simple.index[i]] = F_ref[2]
      uncert_df_simple["F_G_min"][uncert_df_simple.index[i]] = F_G_min[index_uncert]
      uncert_df_simple["F_G_max"][uncert_df_simple.index[i]] = F_G_max[index_uncert]
      uncert_df_simple["F_m_min"][uncert_df_simple.index[i]] = F_m_min[index_uncert]
      uncert_df_simple["F_m_max"][uncert_df_simple.index[i]] = F_m_max[index_uncert]
      uncert_df_simple["F_M_min"][uncert_df_simple.index[i]] = F_M_min[index_uncert]
      uncert_df_simple["F_M_max"][uncert_df_simple.index[i]] = F_M_max[index_uncert]
      uncert_df_simple["F_r_min"][uncert_df_simple.index[i]] = F_r_min[index_uncert]
      uncert_df_simple["F_r_max"][uncert_df_simple.index[i]] = F_r_max[index_uncert]
      uncert_df_simple["F_phi_min"][uncert_df_simple.index[i]] = F_phi_min[index_uncert]
      uncert_df_simple["F_phi_max"][uncert_df_simple.index[i]] = F_phi_max[index_uncert]
      uncert_df_simple["F_z_min"][uncert_df_simple.index[i]] = F_z_min[index_uncert]
      uncert_df_simple["F_z_max"][uncert_df_simple.index[i]] = F_z_max[index_uncert]
      uncert_df_simple["F_R_prime_min"][uncert_df_simple.index[i]] = F_R_prime_min[index_uncert]
      uncert_df_simple["F_R_prime_max"][uncert_df_simple.index[i]] = F_R_prime_max[index_uncert]
   return uncert_df_simple

def uncert_ncal_O4_file_gen(rotors_df):
   G = rotors_df["G"].to_numpy()
   rho = rotors_df["rho"].to_numpy()
   b_rot = rotors_df["b_rot"].to_numpy()
   R_rot = rotors_df["R_rot"].to_numpy()
   r = rotors_df["r"].to_numpy()
   phi = rotors_df["phi"].to_numpy()
   z = rotors_df["z"].to_numpy()
   dG = rotors_df["dG"].to_numpy()
   drho = rotors_df["drho"].to_numpy()
   db_rot = rotors_df["db_rot"].to_numpy()
   dR_rot = rotors_df["dR_rot"].to_numpy()
   dr = rotors_df["dr"].to_numpy()
   dphi = rotors_df["dphi"].to_numpy()
   dz = rotors_df["dz"].to_numpy()
   # ref
   for i in range(len(rotors_df)):
      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i])
      name = "ref_" + name + str(i)
      file_gen(name, data)
      F_ref = force_extract(FROMAGE(directory, name))

      # G
      name, data = ncal_O4_model(G[i]-dG[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i])
      name = "G_min_" + name + str(i)
      file_gen(name, data)
      F_G_min = force_extract(FROMAGE(directory, name))

      name, data = ncal_O4_model(G[i]+dG[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i])
      name = "G_max_" + name + str(i)
      file_gen(name, data)
      F_G_max = force_extract(FROMAGE(directory, name))

      # R_mirror
      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror-dR_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i])
      name = "R_mirror_min_" + name + str(i)
      file_gen(name, data)
      F_R_mirror_min = force_extract(FROMAGE(directory, name))

      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror+dR_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i])
      name = "R_mirror_max_" + name + str(i)
      file_gen(name, data)
      F_R_mirror_max = force_extract(FROMAGE(directory, name))

      # b_mirror 
      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror-db_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i])
      name = "b_mirror_min_" + name + str(i)
      file_gen(name, data)
      F_b_mirror_min = force_extract(FROMAGE(directory, name))

      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror+db_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i])
      name = "b_mirror_max_" + name + str(i)
      file_gen(name, data)
      F_b_mirror_max = force_extract(FROMAGE(directory, name))

      # rho
      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i]-drho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i])
      name = "rho_min_" + name + str(i)
      file_gen(name, data)
      F_rho_min = force_extract(FROMAGE(directory, name))

      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i]+drho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i])
      name = "rho_max_" + name + str(i)
      file_gen(name, data)
      F_rho_max = force_extract(FROMAGE(directory, name))

      # R_rot
      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i]-dR_rot[i], r[i], phi[i], z[i])
      name = "R_rot_min_" + name + str(i)
      file_gen(name, data)
      F_R_rot_min = force_extract(FROMAGE(directory, name))

      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i]+dR_rot[i], r[i], phi[i], z[i])
      name = "R_rot_max_" + name + str(i)
      file_gen(name, data)
      F_R_rot_max = force_extract(FROMAGE(directory, name))

      # b_rot
      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i]-db_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i])
      name = "b_rot_min_" + name + str(i)
      file_gen(name, data)
      F_b_rot_min = force_extract(FROMAGE(directory, name))

      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i]+db_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i])
      name = "b_rot_max_" + name + str(i)
      file_gen(name, data)
      F_b_rot_max = force_extract(FROMAGE(directory, name))

      # r
      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i]-dr[i], phi[i], z[i])
      name = "r_min_" + name + str(i)
      file_gen(name, data)
      F_r_min = force_extract(FROMAGE(directory, name))

      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i]+dr[i], phi[i], z[i])
      name = "r_max_" + name + str(i)
      file_gen(name, data)
      F_r_max = force_extract(FROMAGE(directory, name))

      # phi
      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i]-dphi[i], z[i])
      name = "phi_min_" + name + str(i)
      file_gen(name, data)
      F_phi_min = force_extract(FROMAGE(directory, name))

      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i]+dphi[i], z[i])
      name = "phi_max_" + name + str(i)
      file_gen(name, data)
      F_phi_max = force_extract(FROMAGE(directory, name))

      # z
      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i]-dz[i])
      name = "z_min_" + name + str(i)
      file_gen(name, data)
      F_z_min = force_extract(FROMAGE(directory, name))

      name, data = ncal_O4_model(G[i], rho_mirror, R_mirror, b_mirror, rho[i], b_rot[i], R_inn, R_rot[i], r[i], phi[i], z[i]+dz[i])
      name = "z_max_" + name + str(i)
      file_gen(name, data)
      F_z_max = force_extract(FROMAGE(directory, name))
      

      uncert_df_ncal_O4["F0"][uncert_df_ncal_O4.index[i]] = F_ref[0]
      uncert_df_ncal_O4["F2"][uncert_df_ncal_O4.index[i]] = F_ref[1]
      uncert_df_ncal_O4["A"][uncert_df_ncal_O4.index[i]] = F_ref[2]
      uncert_df_ncal_O4["F_G_min"][uncert_df_ncal_O4.index[i]] = F_G_min[index_uncert]
      uncert_df_ncal_O4["F_G_max"][uncert_df_ncal_O4.index[i]] = F_G_max[index_uncert]
      uncert_df_ncal_O4["F_R_mirror_min"][uncert_df_ncal_O4.index[i]] = F_R_mirror_min[index_uncert]
      uncert_df_ncal_O4["F_R_mirror_max"][uncert_df_ncal_O4.index[i]] = F_R_mirror_max[index_uncert]
      uncert_df_ncal_O4["F_b_mirror_min"][uncert_df_ncal_O4.index[i]] = F_b_mirror_min[index_uncert]
      uncert_df_ncal_O4["F_b_mirror_max"][uncert_df_ncal_O4.index[i]] = F_b_mirror_max[index_uncert]
      uncert_df_ncal_O4["F_r_min"][uncert_df_ncal_O4.index[i]] = F_r_min[index_uncert]
      uncert_df_ncal_O4["F_r_max"][uncert_df_ncal_O4.index[i]] = F_r_max[index_uncert]
      uncert_df_ncal_O4["F_phi_min"][uncert_df_ncal_O4.index[i]] = F_phi_min[index_uncert]
      uncert_df_ncal_O4["F_phi_max"][uncert_df_ncal_O4.index[i]] = F_phi_max[index_uncert]
      uncert_df_ncal_O4["F_z_min"][uncert_df_ncal_O4.index[i]] = F_z_min[index_uncert]
      uncert_df_ncal_O4["F_z_max"][uncert_df_ncal_O4.index[i]] = F_z_max[index_uncert]
      uncert_df_ncal_O4["F_rho_rot_min"][uncert_df_ncal_O4.index[i]] = F_rho_min[index_uncert]
      uncert_df_ncal_O4["F_rho_rot_max"][uncert_df_ncal_O4.index[i]] = F_rho_max[index_uncert]
      uncert_df_ncal_O4["F_R_rot_min"][uncert_df_ncal_O4.index[i]] = F_R_rot_min[index_uncert]
      uncert_df_ncal_O4["F_R_rot_max"][uncert_df_ncal_O4.index[i]] = F_R_rot_max[index_uncert]
      uncert_df_ncal_O4["F_b_rot_min"][uncert_df_ncal_O4.index[i]] = F_b_rot_min[index_uncert]
      uncert_df_ncal_O4["F_b_rot_max"][uncert_df_ncal_O4.index[i]] = F_b_rot_max[index_uncert]
   return uncert_df_ncal_O4


def uncert_calc_simple(uncert_df):
   uncert_df["dF_G"] = np.abs(uncert_df["F_G_max"] - uncert_df["F_G_min"])/2
   uncert_df["dF_m"] = np.abs(uncert_df["F_m_max"] - uncert_df["F_m_min"])/2
   uncert_df["dF_M"] = np.abs(uncert_df["F_M_max"] - uncert_df["F_M_min"])/2
   uncert_df["dF_r"] = np.abs(uncert_df["F_r_max"] - uncert_df["F_r_min"])/2
   uncert_df["dF_phi"] = np.abs(uncert_df["F_phi_max"] - uncert_df["F_phi_min"])/2
   uncert_df["dF_z"] = np.abs(uncert_df["F_z_max"] - uncert_df["F_z_min"])/2
   uncert_df["dF_R_prime"] = np.abs(uncert_df["F_R_prime_max"] - uncert_df["F_R_prime_min"])/2
   
   uncert_df["dtyp_F"] = np.sqrt(uncert_df["dF_G"]**2 + uncert_df["dF_m"]**2 + uncert_df["dF_M"]**2 + uncert_df["dF_r"]**2 + uncert_df["dF_phi"]**2 + uncert_df["dF_z"]**2 + uncert_df["dF_R_prime"]**2)
   uncert_df["dmax_F"] = uncert_df["dF_G"] + uncert_df["dF_m"] + uncert_df["dF_M"] + uncert_df["dF_r"] + uncert_df["dF_phi"] + uncert_df["dF_z"] + uncert_df["dF_R_prime"]
   uncert_df["Etyp"] = uncert_df["dtyp_F"]/uncert_df[uncert_df.columns[index_uncert]]
   uncert_df["Emax"] = uncert_df["dmax_F"]/uncert_df[uncert_df.columns[index_uncert]]
   return uncert_df


def uncert_calc_ncal_O4(uncert_df):
   uncert_df["dF_G"] = np.abs(uncert_df["F_G_max"] - uncert_df["F_G_min"])/2
   uncert_df["dF_R_mirror"] = np.abs(uncert_df["F_R_mirror_max"] - uncert_df["F_R_mirror_min"])/2
   uncert_df["dF_b_mirror"] = np.abs(uncert_df["F_b_mirror_max"] - uncert_df["F_b_mirror_min"])/2
   uncert_df["dF_r"] = np.abs(uncert_df["F_r_max"] - uncert_df["F_r_min"])/2
   uncert_df["dF_phi"] = np.abs(uncert_df["F_phi_max"] - uncert_df["F_phi_min"])/2
   uncert_df["dF_z"] = np.abs(uncert_df["F_z_max"] - uncert_df["F_z_min"])/2
   uncert_df["dF_rho_rot"] = np.abs(uncert_df["F_rho_rot_max"] - uncert_df["F_rho_rot_min"])/2
   uncert_df["dF_R_rot"] = np.abs(uncert_df["F_R_rot_max"] - uncert_df["F_R_rot_min"])/2
   uncert_df["dF_b_rot"] = np.abs(uncert_df["F_b_rot_max"] - uncert_df["F_b_rot_min"])/2

   uncert_df["dtyp_F"] = np.sqrt(uncert_df["dF_G"]**2 + uncert_df["dF_R_mirror"]**2 + uncert_df["dF_b_mirror"]**2 + uncert_df["dF_r"]**2 + uncert_df["dF_phi"]**2 + uncert_df["dF_z"]**2 + uncert_df["dF_rho_rot"]**2 + uncert_df["dF_R_rot"]**2 + uncert_df["dF_b_rot"]**2)
   uncert_df["dmax_F"] = uncert_df["dF_G"] + uncert_df["dF_R_mirror"] + uncert_df["dF_b_mirror"] + uncert_df["dF_r"] + uncert_df["dF_phi"] + uncert_df["dF_z"] + uncert_df["dF_rho_rot"] + uncert_df["dF_R_rot"] + uncert_df["dF_b_rot"]
   uncert_df["Etyp"] = uncert_df["dtyp_F"]/uncert_df[uncert_df.columns[index_uncert]]
   uncert_df["Emax"] = uncert_df["dmax_F"]/uncert_df[uncert_df.columns[index_uncert]]
   return uncert_df

# Analytical model

def Force_max(X, Y):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   F_max = 2*G*m*M * (d_0**2 + R_prime**2)/(d_0**2-R_prime**2)**2 * np.cos(phi) * np.cos(alpha)
   return F_max

def Force_min(X, Y):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   F_min = 2*G*m*M * 1/(d_0**2+R_prime**2) * np.cos(phi) * np.cos(alpha) * np.cos(xi)
   return F_min

def dForce_max_G(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   return Force_max(X, Y)/G*dG

def dForce_max_m(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   return Force_max(X, Y)/m*dm


def dForce_max_M(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   return Force_max(X, Y)/M*dM


def dForce_max_phi(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   return Force_max(X, Y)/np.cos(phi)*np.sin(phi)*dphi

def dForce_max_r(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   dF_max_r = 2*G*m*M*np.cos(phi)\
      *(\
         - 2*r * (d_0**2+3*R_prime**2)/(d_0**2-R_prime**2)**3*np.cos(alpha)\
         + (d_0**2+R_prime**2)/(d_0**2-R_prime**2)**2*z**2/r**3*(np.cos(alpha))**3\
      )*dr
   return dF_max_r

def dForce_max_z(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   dF_max_z = 2*G*m*M*np.cos(phi)\
      *(\
         - 2*z * (d_0**2+3*R_prime**2)/(d_0**2-R_prime**2)**3*np.cos(alpha)\
         + (d_0**2+R_prime**2)/(d_0**2-R_prime**2)**2*z/r**2*(np.cos(alpha))**3\
      )*dz
   return dF_max_z

def dForce_max_R_prime(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   dF_max_R_prime = 2*G*m*M*np.cos(phi)*np.cos(alpha)\
      *(\
         + 2*R_prime * (3*d_0**2+R_prime**2)/(d_0**2-R_prime**2)**3\
      )*dR_prime
   return dF_max_R_prime

def dForce_min_G(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   return Force_min(X, Y)/G*dG

def dForce_min_m(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   return Force_min(X, Y)/m*dm


def dForce_min_M(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   return Force_min(X, Y)/M*dM


def dForce_min_phi(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   return Force_min(X, Y)/np.cos(phi)*np.sin(phi)*dphi

def dForce_min_r(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   dF_max_r = 2*G*m*M*np.cos(phi)\
      *(\
         - 2*r*np.cos(alpha)*np.cos(xi)/(d_0**2+R_prime**2)**2\
         + z**2/r**3*(np.cos(alpha))**3*np.cos(xi)/(d_0**2+R_prime**2)\
         + np.cos(alpha)/(r**2+R_prime**2)*r/np.cos(xi)*(R_prime/(d_0**2+R_prime**2))**2\
      )*dr
   return dF_max_r

def dForce_min_z(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   dF_max_z = 2*G*m*M*np.cos(phi)\
      *(\
         - 2*z*np.cos(alpha)*np.cos(xi)/(d_0**2+R_prime**2)**2\
         + z/r**2*(np.cos(alpha))**3*np.cos(xi)/(d_0**2+R_prime**2)\
         + np.cos(alpha)/(d_0**2+R_prime**2)*z/np.cos(xi)*(R_prime/(d_0**2+R_prime**2))**2\
      )*dz
   return dF_max_z

def dForce_min_R_prime(X, Y, dX):
   (r, phi, z, G, m, M, R_prime) = X ; (d_0, alpha, xi) = Y
   (dr, dphi, dz, dG, dm, dM, dR_prime) = dX
   dF_max_R_prime = 2*G*m*M*np.cos(phi)*np.cos(alpha)\
      *(\
         - 2*R_prime*np.cos(xi)/(d_0**2+R_prime**2)**2 \
         - 1/(d_0**2+R_prime**2)*R_prime/d_0**2*(np.cos(xi))**3\
      )*dR_prime
   return dF_max_R_prime

def dForce_max_max(X, Y, dX):
   return np.abs(dForce_max_G(X, Y, dX)) + \
      np.abs(dForce_max_m(X, Y, dX)) + \
      np.abs(dForce_max_M(X, Y, dX)) + \
      np.abs(dForce_max_phi(X, Y, dX)) + \
      np.abs(dForce_max_r(X, Y, dX)) + \
      np.abs(dForce_max_z(X, Y, dX)) + \
      np.abs(dForce_max_R_prime(X, Y, dX))

def dForce_min_max(X, Y, dX):
   return np.abs(dForce_min_G(X, Y, dX)) + \
      np.abs(dForce_min_m(X, Y, dX)) + \
      np.abs(dForce_min_M(X, Y, dX)) + \
      np.abs(dForce_min_phi(X, Y, dX)) + \
      np.abs(dForce_min_r(X, Y, dX)) + \
      np.abs(dForce_min_z(X, Y, dX)) + \
      np.abs(dForce_min_R_prime(X, Y, dX))

def uncert_calc_analytic(rotors_df, uncert_df):
   for i in range(len(rotors_df)):
      X = (rotors_df["r"][rotors_df.index[i]],\
         rotors_df["phi"][rotors_df.index[i]]*np.pi/180,\
         rotors_df["z"][rotors_df.index[i]],\
         rotors_df["G"][rotors_df.index[i]],\
         rotors_df["m"][rotors_df.index[i]],\
         rotors_df["M"][rotors_df.index[i]],\
         rotors_df["R_prime"][rotors_df.index[i]])
      dX = (rotors_df["dr"][rotors_df.index[i]],\
         rotors_df["dphi"][rotors_df.index[i]]*np.pi/180,\
         rotors_df["dz"][rotors_df.index[i]],\
         rotors_df["dG"][rotors_df.index[i]],\
         rotors_df["dm"][rotors_df.index[i]],\
         rotors_df["dM"][rotors_df.index[i]],\
         rotors_df["dR_prime"][rotors_df.index[i]])
      d_0 = np.sqrt(rotors_df["r"][rotors_df.index[i]]**2 + rotors_df["z"][rotors_df.index[i]]**2)
      alpha = np.arctan(rotors_df["z"][rotors_df.index[i]]/rotors_df["r"][rotors_df.index[i]])
      xi = np.pi/4 - np.arctan((d_0-rotors_df["R_prime"][rotors_df.index[i]])/(d_0+rotors_df["R_prime"][rotors_df.index[i]]))
      Y = (d_0, alpha, xi) 
      uncert_df["F0"][rotors_df.index[i]] = np.abs(Force_max(X, Y) + Force_min(X, Y))/2
      uncert_df["F2"][rotors_df.index[i]] = np.abs(Force_max(X, Y) - Force_min(X, Y))/2
      uncert_df["A"][rotors_df.index[i]] = np.abs(Force_max(X, Y) - Force_min(X, Y))/2
      uncert_df["dF_G"][rotors_df.index[i]] = np.abs(dForce_max_G(X, Y, dX) - dForce_min_G(X, Y, dX))/2
      uncert_df["dF_m"][rotors_df.index[i]] = np.abs(dForce_max_m(X, Y, dX) - dForce_min_m(X, Y, dX))/2
      uncert_df["dF_M"][rotors_df.index[i]] = np.abs(dForce_max_M(X, Y, dX) - dForce_min_M(X, Y, dX))/2
      uncert_df["dF_r"][rotors_df.index[i]] = np.abs(dForce_max_r(X, Y, dX) - dForce_min_r(X, Y, dX))/2
      uncert_df["dF_phi"][rotors_df.index[i]] = np.abs(dForce_max_phi(X, Y, dX) - dForce_min_phi(X, Y, dX))/2
      uncert_df["dF_z"][rotors_df.index[i]] = np.abs(dForce_max_z(X, Y, dX) - dForce_min_z(X, Y, dX))/2
      uncert_df["dF_R_prime"][rotors_df.index[i]] = np.abs(dForce_max_R_prime(X, Y, dX) - dForce_min_R_prime(X, Y, dX))/2

   uncert_df["dtyp_F"] = np.sqrt(uncert_df["dF_G"]**2 + uncert_df["dF_m"]**2 + uncert_df["dF_M"]**2 + uncert_df["dF_r"]**2 + uncert_df["dF_phi"]**2 + uncert_df["dF_z"]**2 + uncert_df["dF_R_prime"]**2)
   uncert_df["dmax_F"] = uncert_df["dF_G"] + uncert_df["dF_m"] + uncert_df["dF_M"] + uncert_df["dF_r"] + uncert_df["dF_phi"] + uncert_df["dF_z"] + uncert_df["dF_R_prime"]
   uncert_df["Etyp"] = uncert_df["dtyp_F"]/uncert_df[uncert_df.columns[index_uncert]]
   uncert_df["Emax"] = uncert_df["dmax_F"]/uncert_df[uncert_df.columns[index_uncert]]
   return uncert_df

uncert_df_simple = uncert_simple_file_gen(rotors_df)
uncert_df_simple = uncert_calc_simple(uncert_df_simple)
uncert_df_analyt = uncert_calc_analytic(rotors_df, uncert_df_analyt)
uncert_df_ncal_O4 = uncert_ncal_O4_file_gen(rotors_df)
uncert_df_ncal_O4 = uncert_calc_ncal_O4(uncert_df_ncal_O4)

uncert_df_analyt.to_csv(csv_path+"uncert_analyt.csv")
uncert_df_simple.to_csv(csv_path+"uncert_simple.csv")
uncert_df_ncal_O4.to_csv(csv_path+"uncert_ncal_O4.csv")

uncert_df_simple = pd.read_csv(csv_path+"uncert_simple.csv", delimiter=",", index_col="NCal")
uncert_df_ncal_O4 = pd.read_csv(csv_path+"uncert_ncal_O4.csv", delimiter=",", index_col="NCal")

# print("ANALYTIC :")
# print(uncert_df_analyt[["A", "dF_G", "dF_m", "dF_M", "dF_r", "dF_phi", "dF_z", "dF_R_prime", "dtyp_F", "dmax_F", "Etyp", "Emax"]])

# print("SIMPLE :")
# print(uncert_df_simple[["A", "dF_G", "dF_m", "dF_M", "dF_r", "dF_phi", "dF_z", "dF_R_prime", "dtyp_F", "dmax_F", "Etyp", "Emax"]])

# print("NCAL O4 :")
# print(uncert_df_ncal_O4[["A", "dF_G", "dF_R_mirror", "dF_b_mirror", "dF_rho_rot", "dF_R_rot", "dF_b_rot", "dF_r", "dF_phi", "dF_z", "dtyp_F", "dmax_F", "Etyp", "Emax"]])

print(uncert_df_ncal_O4[["A", "dtyp_F", "dmax_F", "Etyp", "Emax"]])
