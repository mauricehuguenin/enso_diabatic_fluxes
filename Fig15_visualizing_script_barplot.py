# Purpose: Create a horizontal bar chart with data from my idealised El Niño, 
#          idealised La Niña, and interannual events
#          The goal is to have my table 1 visualised instead of numbers only

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                          Maurice Huguenin-Virchaux                          #
#                       m.huguenin-virchaux@unsw.edu.au                       #
#                           04. 12. 2019, 14:21 AEST                          #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# % preamble ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# I will probably have to hard-code the data here

# links
# https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/horizontal_barchart_distribution.html


# -----------------------------------------------------------------------------
index = 'idealised_EN'
dWWV = -2.5 # total change in WWV over the discharge period
# -----------------------------------------------------------------------------

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties

font = FontProperties()
font.set_family('serif')
font.set_name('Times New Roman')
font.set_size(15)

# Read in data
df = pd.read_csv('/home/z5180028/MSC_thesis/access_matlab_scripts/plot_data.csv', index_col=0)
df = df/100 # Convert to percentages

def custom_plot(example=True):

    def plot_rect(bottom, left, width, height, color = 'C0'):
        ax.add_patch(patches.Rectangle(
                (left, bottom), width, height, linewidth=1, edgecolor=color, facecolor=color))

    # Create figure and axes
    fig,ax = plt.subplots(figsize=(16,6),tight_layout=True,facecolor='w',edgecolor='k')

#    fig, ax = plt.subplots(1)

    # Define axis ticks ticks
    plt.xticks(np.arange(-6,6,1), np.arange(-6,6,1),color=[.192,.211,.584],
               fontsize=15)
    plt.yticks(np.arange(0,5,0.2), np.arange(0,5,0.2),color=[.192,.211,.584],
                fontsize=15)

    ax.set_yticklabels(['Testing','2015/16 event','1997/98 event','1982/83 event','idealised El Niño'],
                       color=[.192,.211,.584])
    # Define axis limits
    plt.ylim(0.05,0.95)
    plt.xlim(-6,6)

    # Move gridlines to the back of the plot
    plt.gca().set_axisbelow(True)

    # Change color of plot edges
    ax.spines['left'].set_color('lightgray')
    ax.spines['right'].set_color('lightgray')
    ax.spines['top'].set_color('lightgray')

    # Hide y axis ticks
#    plt.gca().tick_params(axis='y', colors='w')

    # Turn on gridlines and set color
#    plt.grid(b=True, axis='both', color='lightgray', alpha=0.5, linewidth=1.5)
    # Add lines
    plt.axvline(x=0, c=[.83,.83,.83])
#    plt.axhline(y=0.7, c=[.83,.83,.83])

    # Add x label
    plt.xlabel('Contribution [10$^{14}$ m$^{3}$]', fontsize=15,fontproperties=font,color=[.192,.211,.584])
    # Define color scheme from negative to positive
    colors = [([0.909, 0.301, 0.196]),  # meridional transport
              ([0.192, 0.211, 0.584]),  # Surface forcing
              ([0.337 , 0.690, 1.]),    # Vertical mixing
              ([0., 0.407, 0.215]),     # ITF
              ([0.992, 0.756, 0.435]),  # Surface volume          
              ([0.733, 0., 0.733])]     # numerical mixing 
#    colors = ['firebrick', 'sandybrown', 'navajowhite', 
#              'khaki', 'lightcyan', 'skyblue', 'steelblue']

    # Process data to plot
    try:
        array = [df.iloc[0,:].values, 
                 df.iloc[1,:].values, 
                 df.iloc[2,:].values, 
                 df.iloc[3,:].values]
    except:
        print('Plotting example data')
        example = True

    if example == True:
        # Example data
        array = [np.array([0.2, 0.1, 0.2, 0.2, 0.1, 0.1]),
                 np.array([0.2, 0.1, 0.2, 0.2, 0.1, 0.1]),
                 np.array([0.2, 0.1, 0.2, 0.2, 0.1, 0.1]),
                 np.array([-1.8,-.2,-1.1,+.3, +.1, +.1])] # idealised El Niño

        # Example data column names
        df = pd.DataFrame(columns=['Meridional Transport',
                                   'Indonesian Throughflow',
                                   'Surface Volume',
                                   'Surface forcing',
                                   'Vertical mixing',
                                   'Numerical mixing'])

    # Define function to process input data into rectangle left corner indices
    def process_data(array):
        left = np.zeros_like(array)
        left[0] = -np.sum(array[0:3]) - 3.1
        left[1] = -np.sum(array[1:3]) - 3.1
        left[2] = -np.sum(array[2:3]) - 3.1
        left[3] = 0
        left[4] = np.sum(array[4:4]) + .3
        left[5] = np.sum(array[4:5]) + .3       
        width = array
        return left, width
    
    left = {}
    width = {}
    for i in range(4): # 0,1,2,3
        left[i], width[i] = process_data(array[i])

    # Plot boxes
    height = 0.13
    bottom = 0.135
    for i in range(len(array)):
        for j in range(len(array[i])):
            plot_rect(bottom=bottom+i*0.2, left=left[i][j], 
                      width=width[i][j], height=height, color = colors[j])
            plt.plot(dWWV, .89, 'v', color='k', ms=10, markersize=15)    # Plot marker at total 
                                                 # change in Warm Water Volume

    # Plot category labels
#    plt.text(-1.1,0.9,'Unfavorable', style='italic', 
#             horizontalalignment='left', verticalalignment='center')
#    plt.text(0,0.9,'Neutral', style='italic', 
#             horizontalalignment='center', verticalalignment='center')
#    plt.text(1.1,0.9,'Favorable', style='italic', 
#             horizontalalignment='right', verticalalignment='center')

##     Plot percentages
#    for i in range(len(med)):
#        plt.text(-1,0.2*(i+1),'{0:.0%}'.format(lo[i]), 
#                 horizontalalignment='left', verticalalignment='center')
#        plt.text(0,0.2*(i+1),'{0:.0%}'.format(med[i]), 
#                 horizontalalignment='center', verticalalignment='center')
#        plt.text(1,0.2*(i+1),'{0:.0%}'.format(hi[i]), 
#                 horizontalalignment='right', verticalalignment='center')


# ~~~~~~~~ LEGEND
#    # Create legend
#    fig, ax = plt.subplots(1, figsize=(6,2))
#    plt.axis('off')
#    plt.gca().set_aspect('equal', adjustable='box')
#
#    # Plot colored circles
#    legend_left = [-0.9, -0.6, -0.3, 0, 0.30, 0.6]
#    for i in range(len(colors)):
#        plot_rect(bottom=0, left=legend_left[i], width=0.2, height=0.2, color = colors[i])
#
#    # Plot labels 1-5
#    for i in range(0,5,2):
#        plt.text(-0.8+0.3*i, 0.25, df.columns[i].replace(' ', '\n'), 
#                 horizontalalignment='center', verticalalignment='bottom')
#        plt.text(-0.5+0.3*i, -0.05, df.columns[i+1].replace(' ', '\n'), 
#                 horizontalalignment='center', verticalalignment='top')
#
#    # Plot last label
#    plt.text(1, 0.25, df.columns[5].replace(' ', '\n'), 
#             horizontalalignment='center', verticalalignment='bottom')
#
#    # Plot label title
#    plt.text(-1, 0.1, 'Scale', fontsize=14,
#             horizontalalignment='right', verticalalignment='center')
#
#    plt.gca().autoscale(enable=True, axis='both', tight=None)

custom_plot('example')
#custom_plot(df)

# --- saving as 300 dpi .PNG image in specified folder ------------------------
save = '/home/z5180028/MSC_thesis/access_figures/contribution/'
plt.savefig(save  + 'contribution_figure_' + index, dpi=500, facecolor='w', 
            edgecolor='w', orientation='landscape', papertype=None, 
            format=None, transparent=False, bbox_inches=None, 
            pad_inches=0.1, metadata=None)



# %%
