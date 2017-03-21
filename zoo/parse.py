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


def parse_date(date):
    '''
    Parse (influenza) dates (source: NCBI Genomes, i.e. GenBank) into
    nicely searchable format for insertion in zoo database. Four formats
    are recognized:

    - 1976
    - 2011/01
    - 2011/01/27  # ymd
    - ('NON', 'Unknown', '')

    Usage:

    parse_date('2019/08')
    # {'d': '', 'm': 8, 'y': 2019}
    '''
    record = {'y': None, 'm': None, 'd': None}
    try:  # pandas nan is of type float and won't split
        s = date.split('/')
    except AttributeError:
        return record

    for i in zip(range(3), 'y m d'.split(' ')):
        j, tag = i
        try:  # NON, -N/A-, ...
            record[tag] = int(s[j])
        except (IndexError, ValueError):
            continue
    return record
