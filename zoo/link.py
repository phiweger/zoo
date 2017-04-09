from pyfaidx import Fasta


def link_access(link, _id, func=None):
    '''Given an identifier (e.g. UUID), access seq in supplement fasta file.

    func .. a key function to pass to pyfaidx.Fasta(), e.g.
    func = lambda x: x.split('|')[0]

    The access is random, i.e. quick, given that an index exists. If no index
    is found (file.fai) one is created.
    '''
    try:
        d = link['link']
        if d['type'] == 'filepath':
            fa = Fasta(d['fp'], key_function=func)
            return str(fa[_id])
        else:
            print('This type of link currently not supported.')
            return
    except KeyError:
        print('No link schema provided as input.')
        return
