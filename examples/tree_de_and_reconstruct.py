import networkx as nx
import matplotlib.pyplot as plt
from progressbar import ProgressBar, UnknownLength
import warnings
from zoo.tree import parse_tree, get_leaves


# ignore matplotlib warnings, stackoverflow, 33792478
warnings.filterwarnings("ignore")


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

tree = "('D':0,('C':0,('B':0,'A':0):0,('E':0,('F':0,'G':0):0):0):0);"
g = parse_tree(data=tree, g=nx.DiGraph())

pos = nx.spring_layout(g)
node_labels = {node: node for node in g.nodes()}
nx.draw(g, pos, arrows=True, labels=node_labels)
plt.show()

# decompose
d = {}
for n in get_leaves(g):
    # redundant: don't keep "root" or node itself
    d[n] = nx.shortest_path(g, source='root', target=n)[:-1]

ghat = nx.DiGraph()
bar = ProgressBar(max_value=UnknownLength)
counter = 0
for k, v in d.items():
    v.append(k)
    while len(v) > 1:
        child = v.pop()
        ghat.add_node(child)
        parent = v[-1]
        ghat.add_node(parent)
        ghat.add_edge(parent, child)
    counter += 1
    bar.update(counter)

print('\nTesting isomorphism.')
try:
    assert nx.is_isomorphic(g, ghat)  # slow
    print('Trees g and ghat are isomorphic.')
except AssertionError:
    print('Trees g and ghat NOT isomorphic, something went wrong.')





