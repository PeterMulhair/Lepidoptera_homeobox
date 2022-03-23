import toytree
import toyplot
import toyplot.svg
import toyplot.pdf
import numpy as np
import pandas as pd


'''
Plot proportion of all
TE classes using RM summary
output.
'''

tandDup_sp = []
with open('lep_TE_summ.tsv') as f:
    next(f)
    for line in f:
        lines = line.split('\t')
        sp = lines[0]
        shx_count = lines[2]
        if int(shx_count) > 6:
            tandDup_sp.append(sp)
        

spdata = np.genfromtxt('lep_TE_summ.tsv', delimiter='\t')

tree = toytree.tree("trimmed_sp_tree.nwk")

canvas = toyplot.Canvas(width=1000, height=800)
ax0 = canvas.cartesian(bounds=(50, 200, 50, 600), padding=15, ymin=0, ymax=66)
ax1 = canvas.cartesian(bounds=(225, 325, 50, 600), padding=15, ymin=0, ymax=66)
ax2 = canvas.cartesian(bounds=(360, 800, 50, 600), padding=15, ymin=0, ymax=66)

tree.draw(axes=ax0, width=400, height=600, use_edge_lengths=False, node_sizes=[8 if i.is_leaf() else 0 for i in tree.idx_dict.values()], node_colors=[toytree.colors[1] if (i.is_leaf()) & (i.name in tandDup_sp) else toytree.colors[0] for i in tree.idx_dict.values()], node_style={"stroke": "black", "stroke-width": 1}, tip_labels_style={"-toyplot-anchor-shift": "15px"});
ax0.show = False

ax1.bars(spdata[1:, 1], along='y',color="#fec44f");
ax1.show = True
ax1.y.show = False
ax1.x.ticks.show = True
ax1.x.label.text='Genome size'

colourlist = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#969696", "#bdbdbd", "#d9d9d9", "#f0f0f0"]

ax2.bars(spdata[1:, 3:-1], along='y', color=colourlist);
ax2.show = True
ax2.y.show = False
ax2.x.ticks.show = True
ax2.x.label.text='TE content (% of genome)'


canvas.text(25, 25, "(a)", style={"font-size": "18px"});
canvas.text(300, 25, "(b)", style={"font-size": "18px"});
canvas.text(550, 25, "(c)", style={"font-size": "18px"});

toyplot.svg.render(canvas, 'TE_summary_tree.svg')
