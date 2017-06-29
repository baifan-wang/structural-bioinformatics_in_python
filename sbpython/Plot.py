#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys,os
import time

def plot_torsion_wheel(angles, title, filename, abz=True):

    t1 = time.time()
    label = [r'$\alpha$', r'$\beta$', r'$\gamma$', r'$\delta$', r'$\epsilon$', r'$\zeta$', r'$\chi$']
    fig = plt.figure(figsize=(7, 7), dpi=300)
    fig.subplots_adjust(wspace=0.3)     # To adjust gaps and margins of subplots

    ax1 = plt.subplot(111, polar=True)  # create a ploar plot
    ax1.set_title(title, fontsize=22, y = 1.2)  # set the title of the plot
    ax1.set_theta_direction(-1)         # Set the direction in which theta increases. clockwise, -1, anticlockwise, 1
    ax1.set_theta_zero_location('N')    # Set the location of theta's zero. May be one of "N", "NW", "W", "SW", "S", "SE", "E", or "NE".
    ax1.grid(True, lw=1, linestyle='-', color='black') # to change properties of radial and circular grid-line from dashed to solid
    ax1.spines["polar"].set_visible(False)


    hist_idx = np.arange(0, 360, 1, dtype=np.int)  # Initializing histogram array over 360 Deg using bin-size of 1
    bins = 360
    angle_hist = np.zeros((7,360), dtype=np.int)

    # the observed torsion angle range for A, B and Z DNA. taken from
    a_range = [(180,190,260,320),(140,220),(30,80,140,160),(60,100),(160,250),(250,320),(160,250)]
    b_range = [(270,330),(130,200),(20,80),(70,180),(160,270),(150,210,230,300),(200,300)]
    z_range = [(40,100,150,250),(150,250),(20,90,160,210),(80,160),(180,300),(40,100,280,340),(50,90,180,220)]

    
    for x in range(7):
        if abz == True:
            for i,l in zip([a_range,b_range,z_range],[1,2,4]):
                a = 0
                while a <len(i[x]):
                    for y in range(i[x][a],i[x][a+1]+1):
                        angle_hist[x][y-1]+=l  # The value in angle_hist from A,B,Z are set to 1,2,4,respectively
                                               # the value in overlap part are the sum of the overlap.
                    a+=2
        for l in angles[x]:
            if l is not None:
                angle_hist[x][int(l)-1]=10  # The value in angle_hist from torsion are set to 10

    print("Plotting...")
    width = np.pi/180  # Width of bin in radian
    color_dict = {1:'yellow', 2:'blue', 3:'green', 4:'red', 5:'orange', 6:'magenta', 7:'cyan'}
    for i in range(len(angle_hist)):
        for j in range(len(angle_hist[i])):
            if angle_hist[i][j] == 10:
                plt.bar(hist_idx[j]*width, 1, width=width, bottom=i+1, color='black', linewidth=0)
            elif angle_hist[i][j] !=10 and angle_hist[i][j] !=0:
                plt.bar(hist_idx[j]*width, 1, width=width+0.01, bottom=i+1, \
                    color=color_dict[angle_hist[i][j]], linewidth=0)
            else:
                pass
    # set the data and range to different pattern by ajusting the width

    # Labeling wheel for each angle type
    ax1.yaxis.set_ticklabels(label, fontsize=15)  # Increase font-size of angle around perimeter
    xticks = ax1.get_xticks()*180/np.pi # Get angle label, which are present at the perimeter of wheel and change it to Degree
    xlabel = []

    # Convert angle to string and add Degree symbol
    for x in xticks:
        xlabel.append('{0}$^o$' .format(int(x)))

    # Change postion of angle label to remove any overlap
    ax1.set_thetagrids(xticks, frac=1.1)
#    theta_angles = np.array([0,90,180,270])
#    ax1.set_thetagrids(theta_angles,visible= True)

    # At last change the fontsize
    ax1.set_rlabel_position(0)  #set the position of the angle label using angle as coordinate
    ax1.xaxis.grid(True, alpha = 0.1)  #turn on the radial grids.
    ax1.xaxis.set_ticklabels(xlabel, fontsize=18)

    # Clean memory
    del angles

    #plt.show()
    fig.savefig(filename, dpi=180, bbox_inches='tight')
    t2 = time.time()
    print('\nTotal time used :%2f s' %(t2-t1))


def plot_phi_psi(torsion, title, filename):
    t1 = time.time()

    phi = [i[0] for i in torsion]
    psi = [i[1] for i in torsion]
    fig = plt.figure(figsize=(7, 7), dpi=300)
    fig.subplots_adjust(wspace=0.3)     # To adjust gaps and margins of subplots

    ax1 = plt.subplot(111)  # create a ploar plot
    ax1.set_title(title, fontsize=22)  # set the title of the plot
    ax1.set_xlim((-180,180))
    ax1.set_ylim((-180,180))
    ax1.set_xlabel('Phi (degrees)',fontsize=14)
    ax1.set_ylabel('Psi (degrees)',fontsize=14)
    major_ticks = [-180,-135,-90,-45,0,45,90,135,180]
    minor_ticks = [0]

    ax1.set_xticks(major_ticks)
    ax1.set_xticks(minor_ticks, minor=True)
    ax1.set_yticks(major_ticks)
    ax1.set_yticks(minor_ticks, minor=True)
    ax1.grid(which='minor', linestyle = '--', linewidth = 1)

    ax1.set_xticklabels(major_ticks, fontsize=12)
    ax1.set_yticklabels(major_ticks, fontsize=12)
    plt.plot(phi, psi, marker='s', linewidth=0, markersize=5)

    fig.savefig(filename, dpi=180, bbox_inches='tight')

    t2 = time.time()
    print('\nTotal time used :%2f s' %(t2-t1))
    pass