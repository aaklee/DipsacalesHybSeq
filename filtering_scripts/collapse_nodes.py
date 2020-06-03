#!/usr/bin/python

import sys
import os
import dendropy

tre_f = os.path.abspath(sys.argv[1])
bootstrap = int(sys.argv[2])

data = tre_f.split('/')[-1].split('_')[0]

tre_list = []
with open(tre_f, 'r') as inf:
    for line in inf:
        tre_list.append(line.strip())


collapsed_trees = []
for f in tre_list:
    collapse = 0
    locus = f.split('/')[-1].split('.')[1]

    tre = dendropy.Tree.get(path=f, schema='newick')
    for edge in tre.postorder_edge_iter():
        if edge.head_node.label and int(edge.head_node.label) < bootstrap:
            edge.collapse()
            collapse += 1

    print('{}, collapsed {} nodes with bootstrap < {}'.format(locus, collapse, bootstrap))

    out = '{}.{}_collapsed.gene.tre'.format(data, locus)
    collapsed_trees.append(out)

    with open(out, 'w+') as outf:
        tre.write(file=outf, schema='newick')

os.system('cat {} > {}_collapsed.gene.tre'.format(' '.join(collapsed_trees), data))
