from scan_track import (ReadTrack, Box)
from scan_track import (stratification, neighbourhood)
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    track = ReadTrack("/home/imc/ksen/stretching/")
    # box = Box(100, 100, 100)
    n_bins = 40
    rho = np.zeros(n_bins)
    rho_h = np.zeros(n_bins)
    rho_t = np.zeros(n_bins)
    rho_w = np.zeros(n_bins)
    counter = 0
    while track.one_step():
        counter += 1
        rho_head = np.zeros(n_bins)
        rho_tail = np.zeros(n_bins)
        rho_water = np.zeros(n_bins)
        box_z = track.box.z
        v = box_z/n_bins*track.box.x*track.box.y
        z = track.z
        z = (z+box_z/2)/box_z*n_bins
        for i, t in enumerate(track.btype):
            n_layer = int(z[i]) 
            if t == 1:
                rho_head[n_layer] += 1
            if t == 2:
                rho_tail[n_layer] += 1
            if t == 3:
                rho_water[n_layer] += 1
        z0 = np.argmax(rho_tail)
        rho_t += np.hstack([rho_tail[z0:], rho_tail[:z0]]) 
        rho_h += np.hstack([rho_head[z0:], rho_head[:z0]]) 
        rho_w += np.hstack([rho_water[z0:], rho_water[:z0]]) 
    rho = rho_t + rho_h + rho_w
    z_layers = np.linspace(0,box_z, n_bins) 
    
    plt.plot(z_layers,rho_h/v/counter,label='head')
    plt.plot(z_layers,rho_t/v/counter, label='tail')
    plt.plot(z_layers,rho_w/v/counter,label='water')
    plt.plot(z_layers,rho/v/counter, label='all')
    plt.legend()
    plt.xlabel('$z$')
    plt.ylabel(r'$\rho(z)$')
    plt.show()    
        
        # n = len(x)
        # energy = 0.0
        # for i in range(n-1):
        #     for j in range(i+1, n):
        #         dx = x[i] - x[j]
        #         dy = y[i] - y[j]
        #         dz = z[i] - z[j]
        #         dx, dy, dz = track.box.periodic_correct(dx, dy, dz)
        #         energy += dx**2+dy**2+dz**2
        # print(energy)       
        
        # print(np.mean(fz**2), np.mean(fx**2+fy**2), np.mean(fx**2))
       
        # print(track.box.z)
        # z_coords = dict()
        # for zz in z:
        #     if int(zz) in z_coords:
        #         z_coords[int(zz)] += 1
        #     else:
        #         z_coords[int(zz)] = 1
        # keys = sorted(z_coords)
        # # print(keys[-1]-keys[0]+1, len(keys))
        # if abs(keys[0]) != keys[-1] or (keys[-1]-keys[0]+1) == len(keys):
        #     print(track.time_step, "hole")
        # # for key in keys:
        # #     print(key, z_coords[key])n = len(x)
        # energy = 0.0
        # for i in range(n-1):
        #     for j in range(i+1, n):
        #         dx = x[i] - x[j]
        #         dy = y[i] - y[j]
        #         dz = z[i] - z[j]
        #         dx, dy, dz = track.box.periodic_correct(dx, dy, dz)
        #         energy += dx**2+dy**2+dz**2
        # print(energy)
        # # clusters = stratification(len(x), neighbourhood(x, y, z, box=box, radius=1.0))
        # # print(len(clusters))
        
        # # break