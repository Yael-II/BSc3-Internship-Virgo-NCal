"""
FETA : Fromage Element Tracker and Animator
Yaël MOUSSOUNI
University of Strasbourg, Department of Physics and Engineering
"""
# Imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# Variables and parameters
directory = "./out/" # directory of output files
name_rotor = "rotor" # rotor position files names (without _N.txt)
name_mirror = "mirror.txt" # mirror position file name (with extension)
name_disp = "disp.txt" # displacement file name (with extension)

printed_text = "" # optional text to add to the figure (leave blank "" for default)

view_infos = True # include informations ? (ax3)
view_text = True # include text ? (ax4)
to_gif = True # export to GIF ?
to_screen = True # show the result ?

# PLT Setup
plt.style.use("default")
fig = plt.figure()
fig.set_size_inches(12.0, 7.0)
ax1 = plt.subplot2grid((3,3), (0,0), projection='3d', rowspan=3, colspan=2) # 3D View
ax2 = plt.subplot2grid((3,3), (1,2), rowspan=1, colspan=1) # graph
ax3 = plt.subplot2grid((3,3), (0,2), rowspan=1, colspan=1) # other informations
ax4 = plt.subplot2grid((3,3), (2,2), rowspan=1, colspan=1) # other informations

ax2.grid()
# text zones (no grid, no axis, nothing but text)
ax3.get_xaxis().set_visible(False)
ax3.get_yaxis().set_visible(False)
ax4.get_xaxis().set_visible(False)
ax4.get_yaxis().set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['bottom'].set_visible(False)
ax4.spines['left'].set_visible(False)
ax3.set_xlim(0,1)
ax3.set_ylim(0,1)
ax4.set_xlim(0,1)
ax4.set_ylim(0,1)
# Data reading
file_mirror = open(directory+name_mirror, "r")
file_mirror.readline() ; file_mirror.readline()
data_mirror = file_mirror.readlines()
file_mirror.close()

file_disp = open(directory+name_disp, "r")
file_disp.readline() ; file_disp.readline()
data_disp = file_disp.readlines()
file_disp.close()

list_file = os.listdir(directory)
list_file = [f for f in list_file if name_rotor in f]
data_rotor = []
for i in range(len(list_file)) :
      file_rotor = open(directory+name_rotor+"_"+str(i)+".txt", "r")
      file_rotor.readline() ; file_rotor.readline()
      data_rotor.append(file_rotor.readlines())
      file_rotor.close()

# Data manipulation
## Mirror
x_mirror = []
y_mirror = []
z_mirror = []
m_mirror = []
for i in range(len(data_mirror)):
   x, y, z, m = data_mirror[i].replace(" \n", "").split(" ")
   x_mirror.append(float(x))
   y_mirror.append(float(y))
   z_mirror.append(float(z))
   m_mirror.append(float(m))
x_mirror = np.array(x_mirror)
y_mirror = np.array(y_mirror)
z_mirror = np.array(z_mirror)
m_mirror = np.array(m_mirror)

## Rotor
x_rotor = []
y_rotor = []
z_rotor = []
m_rotor = []

for i in range(len(data_rotor)):
   x_rotor.append([])
   y_rotor.append([])
   z_rotor.append([])
   m_rotor.append([])
   for j in range(len(data_rotor[i])):
      x, y, z, m = data_rotor[i][j].replace(" \n", "").split(" ")
      x_rotor[i].append(float(x))
      y_rotor[i].append(float(y))
      z_rotor[i].append(float(z))
      m_rotor[i].append(float(m))
x_rotor = np.array(x_rotor)
y_rotor = np.array(y_rotor)
z_rotor = np.array(z_rotor)
m_rotor = np.array(m_rotor)

theta = []
disp = []

for i in range(len(data_disp)):
   t, d =  data_disp[i].replace(" \n", "").split('    ')
   theta.append(float(t))
   disp.append(float(d))
# Maximum masses
m_min = np.min([np.min(m_rotor[0]), np.min(m_mirror)])
m_max = np.max([np.min(m_rotor[0]), np.min(m_mirror)])
m_tot_mirror = np.sum(m_mirror)
m_tot_rotor = np.sum(m_rotor[0])
# Aninmation
def anim(n, scat_rotor, plot_curr):
   plot_curr[0].set_data((theta[n], disp[n]))
   scat_rotor._offsets3d = (x_rotor[n], y_rotor[n], z_rotor[n])
   return scat_rotor, plot_curr
# Plotting
scat_mirror = ax1.scatter(x_mirror, y_mirror, z_mirror, c=m_mirror, cmap='cividis_r', vmin=m_min, vmax=m_max)
scat_rotor = ax1.scatter(x_rotor[0], y_rotor[0], z_rotor[0], c=m_rotor[0], cmap='cividis_r', vmin=m_min, vmax=m_max)
plot_disp = ax2.plot(theta, disp, "+--", color='#3898ff')
plot_curr = ax2.plot(theta[0], disp[0], "o-", color='#003763')
plt.suptitle("FETA, a FROMAGE extension")
ax1.set_title("Position of elements")
ax1.set_xlabel("$x$/m")
ax1.set_ylabel("$y$/m")
ax1.set_zlabel("$z$/m")
ax2.set_title("Mirror displacement with rotation of the rotor")
ax2.set_xlabel("Rotor angle $\\theta$/deg")
ax2.set_ylabel("Mirror displacement $\\delta$/m")
ax1.set_aspect('equal')
# ax1.set_xlim(ax1.get_xlim())
# ax1.set_ylim(ax1.get_ylim())
# ax1.set_zlim(ax1.get_zlim())
# ax2.set_xlim(ax2.get_xlim())
# ax2.set_ylim(ax2.get_ylim())

plt.colorbar(scat_rotor, cmap='cividis_r', orientation='horizontal', label='Mass of the element $m_i$/kg', shrink=0.5)

if view_infos :
   ax3.text(0,0.5, 'Total mirror mass : {} kg\nTotal rotor mass : {} kg'.format(str(np.round(m_tot_mirror, 5)), str(np.round(m_tot_rotor, 5))), horizontalalignment='left', verticalalignment='center', transform=ax3.transAxes)

if view_text and printed_text == "" :
   ax4.text(0.5,0.5, """ "You spin me right 'round, baby, right 'round" """, style='italic', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes)
   ax4.text(1,0.3, """ Dead or alive """, horizontalalignment='right', verticalalignment='center', transform=ax4.transAxes)
elif view_text :
   ax4.text(0.5,0.5, printed_text, horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes)
ani = animation.FuncAnimation(fig, anim, len(theta), fargs=(scat_rotor, plot_curr), interval=100)

plt.tight_layout()
## Output
writer = animation.PillowWriter(fps=10, metadata=dict(artist='FETA'), bitrate=1800)
if to_gif : ani.save('feta.gif', writer=writer)
plt.savefig("feta.png")
if to_screen : plt.show()

