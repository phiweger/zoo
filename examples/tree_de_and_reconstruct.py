import networkx as nx
import matplotlib.pyplot as plt
import re


def get_leaves(g):
    '''Return leaves (terminal nodes) of an DAG, stackoverflow, 31946253.'''
    return [x for x in g.nodes_iter() if g.out_degree(
        x) == 0 and g.in_degree(x) == 1]


tree = '(D,(C,(B,A),(E,(F,G))))'


'''
g
'''


g = nx.DiGraph()
index = 0
queue = []  # FIFO

for i in tree:
    if i == '(':
        if index == 0:
            queue.append('root')
        else:
            g.add_node(str(index))  # is_isomorphic does not work with str/ int
            g.add_edge(queue[-1], str(index))
            queue.append(str(index))
        index += 1
    elif i == ',':
        pass
    elif i == ')':
        _ = queue.pop()
    else:  # letter
        g.add_node(i)  # leaf
        g.add_edge(queue[-1], i)  # directed edge, from -> to
try:
    assert not queue
    print('Tree built.')
except AssertionError:
    print('The queue should be empty, something went wrong.')


# seems right
pos = nx.spring_layout(g)
node_labels = {node: node for node in g.nodes()}
nx.draw(g, pos, arrows=True, labels=node_labels)
plt.show()


n = 'C'
print([i for i in g.predecessors_iter(n)])
print([i for i in g.successors_iter(n)])


nx.shortest_path(g, source='root', target='A')
'''
['root', 1, 2, 'A']
save ['root', 1, 2]
This is robust to removal of entries as well, although the branch weights
will be different in this case I think. Maybe issue warning and that's that.
'''


'''
ghat
'''


d = {}
for n in get_leaves(g):
    # redundant: don't keep "root" or node itself
    d[n] = nx.shortest_path(g, source='root', target=n)[:-1]

ghat = nx.DiGraph()
for k, v in d.items():
    v.append(k)
    while len(v) > 1:
        child = v.pop()
        ghat.add_node(child)
        parent = v[-1]
        ghat.add_node(parent)
        ghat.add_edge(parent, child)
try:
    assert nx.is_isomorphic(g, ghat)
    print('Trees g and ghat isomorphic.')
except AssertionError:
    print('Trees g and ghat NOT isomorphic, something went wrong.')


# --------------------------------------------------------------------------
# below is under development


def parse_tree(treestr, tree, verbose=False):
    '''
    tree = "('D':0,('C':0,('B':0,'A':0):0,('E':0,('F':0,'G':0):0):0):0)"

    for some guidance, see:
    https://github.com/blab/baltic/blob/master/baltic.py

    treestr .. a tree string
    tree .. a networkx graph in the shape of a tree

    https://pyformat.info/
    '''

    tree = nx.DiGraph()
    i = 0
    stored_i = None

    while i < len(treestr):
        if stored_i == i and verbose is True:
            print('{} >{}<'.format(i, treestr[i]))

        print(i, stored_i)

        assert (stored_i != i), '\nTree string unparseable.\nStopped at >>{}<<\nstring region looks like this: {}'.format(treestr[i], treestr[i:i+5000]) ## make sure that you've actually parsed something last time, if not - there's something unexpected in the tree string
        stored_i = i

        if treestr[i] == '(':
            if verbose is True:
                print('adding node {}'.format(i))
            tree.add_node(i)
            i += 1

        cerberus = re.match('(\(|,)([0-9]+)(\[|\:)', treestr[i-1:i+100])
        print(cerberus)
        if cerberus is not None:
            if verbose is True:
                print('{} adding leaf (BEAST) {}'.format(i, cerberus.group(2)))
            tree.add_leaf(i, cerberus.group(2)) ## add tip
            i += len(cerberus.group(2)) ## advance in tree string by however many characters the tip is encoded


# rewrite baltic.make_tree in Python 3 and w/ networkx graph structure.














