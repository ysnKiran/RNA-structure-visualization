import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D

f = open('data.txt')

color_codes = {
    'A': 'red',
    'C': 'blue',
    'G': 'orange',
    'U': 'green'
}

seq = f.readline()
seq = seq.strip('\n')
len = len(seq.strip('\n'))
fig, ax = plt.subplots()
for i in range(len):
    if i:
        line = patches.Arc((i - 0.5, 0), 0.75, 0.01, 0, 180, 360)
        ax.add_patch(line)
    plt.scatter(i, 0, c=color_codes[seq[i]])
    plt.annotate(seq[i], (i, 0), (i - (len * 0.005), -len * 0.1))

pairs = int(f.readline())
for i in range(pairs):
    idx1, idx2, base1, base2 = f.readline().split(' ')
    idx1 = int(idx1)
    idx2 = int(idx2)
    semicircle = patches.Arc(((idx1 + idx2)/2, 0), idx2 - idx1, (idx2 - idx1)/2, 0, 0, 180, fill=False, color='red', linestyle='--')
    ax.add_patch(semicircle)
    plt.scatter(idx1, 0, c=color_codes[base1.strip('\n')], edgecolors='black')
    plt.scatter(idx2, 0, c=color_codes[base2.strip('\n')], edgecolors='black')

plt.axis('off')
plt.title("Secondary Structure for RNA Sequence")
plt.xlim(-1, len)
plt.ylim(-len/2, len/2)
plt.show()

fig_name = ('.\\visualization\\examples\\' + seq + '.png').strip('\n')
fig.savefig(fig_name, dpi=600)