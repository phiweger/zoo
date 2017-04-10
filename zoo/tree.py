import networkx as nx
import re


def get_leaves(g):
    '''Return leaves (terminal nodes) of an DAG, stackoverflow, 31946253.'''
    return [x for x in g.nodes_iter() if g.out_degree(
        x) == 0 and g.in_degree(x) == 1]


def parse_tree(data, g=nx.DiGraph(), verbose=False):
    '''
    tree .. a tree string
    g .. a networkx graph in the shape of a tree

    tree = "('D':0,('C':0,('B':0,'A':0):0,('E':0,('F':0,'G':0):0):0):0);"
    g = nx.DiGraph()

    We need a directed graph here, e.g. because we identify leaves as
    nodes w/ in-degree 1 and out-degree 0.

    For some guidance on the tree parsing, see the excellent:
    https://github.com/blab/baltic/blob/master/baltic.py

    For information in legacy Python2 format (for rewrite of baltic):
    https://pyformat.info/
    '''
    i = 0
    # is an adjustable index along the tree string, it is incremented to
    # advance through the string
    stored_i = None
    # store the i at the end of the loop, to make sure we haven't gotten stuck
    # somewhere in an infinite loop

    index = 0
    queue = []  # FIFO

    while i < len(data):
        # while there's characters left in the tree string - loop away
        if stored_i == i and verbose is True:
            print('{} >{}<'.format(i, data[i]))

        assert (stored_i != i), '''\nTree string unparseable\nStopped at \
        >>{}<<\nstring region looks like this: {}'''.format(
            data[i], data[i:i+5000])
        # make sure that you've actually parsed something last time, if not -
        # there's something unexpected in the tree string
        stored_i = i  # store i for later

        '''
        case 1
        '''
        if data[i] == '(':  # look for new nodes
            if verbose is True:
                print('{} adding node'.format(i))
            # ll.add_node(i)  # add node to current node in tree ll
            i += 1  # advance in tree string by one character
            if index == 0:
                queue.append('root')
            else:
                g.add_node(str(index))
                # is_isomorphic does not work with str/ int
                g.add_edge(queue[-1], str(index))
                queue.append(str(index))
            index += 1

        '''
        case 2
        '''
        cerberus = re.match(
            '(\(|,)(\'|\")*([A-Za-z\_\-\|\.0-9\?\/]+)(\'|\"|)(\[)*',
            data[i-1:i+200]
            )
        # look for tips with unencoded names - if the tips have some
        # unusual format you'll have to modify this
        if cerberus is not None:
            if verbose is True:
                print('{} adding leaf (non-BEAST) {}'.format(
                    i, cerberus.group(3)))

            i += len(cerberus.group(3)) + \
                cerberus.group().count("'") + cerberus.group().count('"')
            # advance in tree string by however many characters the tip is
            # encoded

            g.add_node(cerberus.group(3))  # leaf
            g.add_edge(queue[-1], cerberus.group(3))
            # directed edge, from -> to

        '''
        case 3
        '''
        microcerberus = re.match('(\:)*([0-9\.\-Ee]+)', data[i:i+100])
        # look for branch lengths without comments
        if microcerberus is not None:
            if verbose is True:
                print('adding branch length ({}) {}'.format(
                    i, float(microcerberus.group(2))))
            # length of current node
            i += len(microcerberus.group())
            # advance in tree string by however many characters it took to
            # encode branch length

        '''
        case 4
        '''
        if data[i] == ')':
            queue.pop()

        if data[i] == ',' or data[i] == ')':
            # look for bifurcations or clade ends
            i += 1  # advance in tree string
            # ll.cur_node=ll.cur_node.parent

        if data[i] == ';':  # look for string end
            return g  # end loop


'''
def make_tree(data,ll,verbose=False):
    """
    data is a tree string, ll (LL) is an instance of a tree object
    """
    i=0 ## is an adjustable index along the tree string, it is incremented to advance through the string
    stored_i=None ## store the i at the end of the loop, to make sure we haven't gotten stuck somewhere in an infinite loop
    
    while i < len(data): ## while there's characters left in the tree string - loop away
        if stored_i == i and verbose==True:
            print('{} >{}<'.format(i,data[i]))
        
        assert (stored_i != i),'\nTree string unparseable\nStopped at >>{}<<\nstring region looks like this: {}'.format(data[i],data[i:i+5000]) ## make sure that you've actually parsed something last time, if not - there's something unexpected in the tree string
        stored_i=i ## store i for later
        
        if data[i] == '(': ## look for new nodes
            if verbose==True:
                print('{} adding node'.format(i))
            ll.add_node(i) ## add node to current node in tree ll
            i+=1 ## advance in tree string by one character
            
        cerberus=re.match('(\(|,)([0-9]+)(\[|\:)',data[i-1:i+100]) ## look for tips in BEAST format (integers).
        if cerberus is not None:
            if verbose==True:
                print('{} adding leaf (BEAST) {}'.format(i,cerberus.group(2)))
            ll.add_leaf(i,cerberus.group(2)) ## add tip
            i+=len(cerberus.group(2)) ## advance in tree string by however many characters the tip is encoded
            
        cerberus=re.match('(\(|,)(\'|\")*([A-Za-z\_\-\|\.0-9\?\/]+)(\'|\"|)(\[)*',data[i-1:i+200])  ## look for tips with unencoded names - if the tips have some unusual format you'll have to modify this
        if cerberus is not None:
            if verbose==True:
                print('{} adding leaf (non-BEAST) {}'.format(i,cerberus.group(3)))
            ll.add_leaf(i,cerberus.group(3).strip('"').strip("'"))  ## add tip
            i+=len(cerberus.group(3))+cerberus.group().count("'")+cerberus.group().count('"') ## advance in tree string by however many characters the tip is encoded

        cerberus=re.match('\)([0-9]+)\[',data[i-1:i+100]) ## look for multitype tree singletons.
        if cerberus is not None:
            if verbose==True:
                print('{} adding multitype node {}'.format(i,cerberus.group(1)))
            i+=len(cerberus.group(1))

        cerberus=re.match('(\:)*\[&([A-Za-z\_\-{}\,0-9\.\%=\"\+!#]+)\]',data[i:])## look for MCC comments
        #cerberus=re.match('\[&[A-Za-z\_\-{}\,0-9\.\%=\"\+]+\]',data[i:])## look for MCC comments
        if cerberus is not None:
            if verbose==True:
                print('{} comment: {}'.format(i,cerberus.group(2)))
            comment=cerberus.group(2)
            numerics=re.findall('[A-Za-z\_\.0-9]+=[0-9\-Ee\.]+',comment) ## find all entries that have values as floats
            strings=re.findall('[A-Za-z\_\.0-9]+=["|\']*[A-Za-z\_0-9\.\+]+["|\']*',comment) ## strings
            treelist=re.findall('[A-Za-z\_\.0-9]+={[A-Za-z\_,{}0-9\.]+}',comment) ## complete history logged robust counting (MCMC trees)
            sets=re.findall('[A-Za-z\_\.0-9\%]+={[A-Za-z\.\-0-9eE,\"\_]+}',comment) ## sets and ranges
            figtree=re.findall('\![A-Za-z]+=[A-Za-z0-9#]+',comment)
            
            for vals in strings:
                tr,val=vals.split('=')
                if '+' in val:
                    val=val.split('+')[0] ## DO NOT ALLOW EQUIPROBABLE DOUBLE ANNOTATIONS (which are in format "A+B") - just get the first one
                ll.cur_node.traits[tr]=val.strip('"')
            
            for vals in numerics: ## assign all parsed annotations to traits of current branch
                tr,val=vals.split('=') ## split each value by =, left side is name, right side is value
                ll.cur_node.traits[tr]=float(val)
                
            for val in treelist:
                tr,val=val.split('=')
                microcerberus=re.findall('{([0-9]+,[0-9\.\-e]+,[A-Z]+,[A-Z]+)}',val)
                ll.cur_node.traits[tr]=[]
                for val in microcerberus:
                    codon,timing,start,end=val.split(',')
                    ll.cur_node.traits[tr].append((int(codon),float(timing),start,end))
                
            for vals in sets:
                tr,val=vals.split('=')
                if 'set' in tr:
                    ll.cur_node.traits[tr]=[]
                    for v in val[1:-1].split(','):
                        if 'set.prob' in tr:
                            ll.cur_node.traits[tr].append(float(v))
                        else:
                            ll.cur_node.traits[tr].append(v.strip('"'))
                elif 'range' in tr or 'HPD' in tr:
                    ll.cur_node.traits[tr]=map(float,val[1:-1].split(','))
                else:
                    print('some other trait: {}'.format(vals))
            
            if len(figtree)>0:
                print('FigTree comment found, ignoring')
            
            i+=len(cerberus.group()) ## advance in tree string by however many characters it took to encode labels
            
        cerberus=re.match('([A-Za-z\_\-0-9\.]+)(\:|\;)',data[i:])## look for old school node labels
        if cerberus is not None:
            if verbose==True:
                print('old school comment found: {}'.format(cerberus.group(1)))
            ll.cur_node.traits['label']=cerberus.group(1)
            
            i+=len(cerberus.group(1))
            
        microcerberus=re.match('(\:)*([0-9\.\-Ee]+)',data[i:i+100]) ## look for branch lengths without comments
        if microcerberus is not None:
            if verbose==True:
                print('adding branch length ({}) {}'.format(i,float(microcerberus.group(2))))
            ll.cur_node.length=float(microcerberus.group(2)) ## set branch length of current node
            i+=len(microcerberus.group()) ## advance in tree string by however many characters it took to encode branch length

        if data[i] == ',' or data[i] == ')': ## look for bifurcations or clade ends
            i+=1 ## advance in tree string
            ll.cur_node=ll.cur_node.parent

        if data[i] == ';': ## look for string end
            break ## end loop
'''
