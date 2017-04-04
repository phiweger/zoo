import re


def parse_nomenclature_iav(notation):
    '''
    parse the data from influenza name into dict

    test:
    a = 'A/mallard/Interior Alaska/6MP0155/2006(H3N8)'
    b = 'A/blue-winged Teal/Minnesota/AI09-2977/2009(H4N8)'
    c = 'A/Peru/PER175/2012(H3N2)'
    '''

    # regexp structure check

    feat = notation.split('/')
    assert feat[0] == 'A', 'Not an influenza virus of type "A".'

    try:
        subtype = re.search(r'\((.*?)\)', notation).group(1)
        feat[-1] = feat[-1].split('(')[0]
        if len(feat) == 4:
            host = 'human'
            geo = feat[1]
        else:
            host = feat[1].lower()
            geo = feat[2]

        return {
            'type': feat[0],
            'host': host,
            'subtype': subtype,
            'location': geo
            }
    except IndexError:
        return {
            'type': None,
            'host': None,
            'subtype': None,
            'location': None
            }


def construct_nomenclature_iav(d):
    '''Transform dict to IAV nomenclature.'''
    pass
