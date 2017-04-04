import re
from Bio.SeqFeature import FeatureLocation, CompoundLocation


def parse_date(date):
    '''
    Parse (influenza) dates (source: NCBI Genomes, i.e. GenBank) into
    nicely searchable format for insertion in zoo database. Four formats
    are recognized:

    - 1976
    - 2011/01
    - 2011/01/27  # ymd
    - ('NON', 'Unknown', '', ...)

    Note that year, month, day order is assumed.

    Usage:

    parse_date('2019/08')
    # {'d': '', 'm': 8, 'y': 2019}
    '''
    record = {'y': None, 'm': None, 'd': None}
    try:
        s = re.split('-|/', date)
    except TypeError:  # date(3.4) or pandas nan of type float
        return record

    for i in zip(range(3), 'y m d'.split(' ')):
        j, tag = i
        try:  # NON, -N/A-, ...
            record[tag] = int(s[j])
        except (IndexError, ValueError):
            continue
    return record


# DEPRECATED, kept for reference
# def parse_location(location, source='genbank'):
#     '''
#     Example:

#     parse_gb_location('1->1374')
#     # ('1', '1374', 1)  # start, end, fuzzy
#     '''
#     if source == 'genbank':
#         start, end = re.findall('\\d+', location)
#         fuzzy = int(any([i in location for i in ['>', '<']]))
#         return start, end, fuzzy
#     raise AttributeError('Unknown format.')


def location_tostr(loc):
    '''Parse FeatureLocation or dict to str.

    In: An iterable (list, generator) of either

    - Biopython FeatureLocation and/ or CompoundLocation objects
    - dicts

    Out: A generator of location strings in the format '[20:>30](?)'.

    Example dict format (output from zoo.parse.location_todict():

    {
        'start': 10,
        'end': 40,
        'strand': '?',
        'pstart': 0,
        'pend': 1
    }

    Example:

    loc = ['[5:10](-)', '[20:>30](?)']
    gen = location_todict(loc)
    assert next(location_tostr(gen)) == loc[0]  # True
    '''
    for i in loc:
        if isinstance(i, (FeatureLocation, CompoundLocation)):
            if len(i.parts) == 1:
                yield str(i)
            else:
                yield [str(j) for j in i.parts]

        elif isinstance(i, dict):
            try:
                ds = sorted(i.items())
                end, pend, pstart, start, strand = [v[1] for v in ds]
            except ValueError:  # too many or too little keys
                print('Make sure loc is in the right format.')
                return
            pstart = '<' if pstart else ''
            pend = '>' if pend else ''

            yield '[{}{}:{}{}]({})'.format(
                pstart, start, pend, end, strand
                )
        else:
            print('Currently only supports FeatureLocation and dict types.')
            return


def location_todict(loc):
    '''Parse location str to dict.

    In: An iterable of (location) strings in format '[20:>30](?)'.
    Out: An iterable of dicts.

    Example:

    loc = ['[5:10](-)', '[20:>30](?)']
    gen = location_todict(loc)
    next(gen)['start'] == 5
    assert {v for k, v in next(gen)} == {20, 30, 1, 0, '?'}

    '''
    for i in loc:
        start, end = re.search(r'\[(.*?):(.*?)\]', i).group(1, 2)
        strand = re.search(r'\((.)\)', i).group(1)
        yield {
            'start': int(start.replace('<', '')),
            'end': int(end.replace('>', '')),
            'strand': strand,
            'pstart': int('<' in start),  # p .. partial
            'pend': int('>' in end)       # 0 no, 1 yes
            }


def location_extract(loc, seq):
    '''Extract subsequence given coordinates in location argument.

    In: A location dict, as generated from zoo.parse.location_todict()
    (see format requirement in corresponding function doc) and a sequence
    string like "ACTAGANNCAT...".
    Out: A generator of subsequence strings.

    We could concatenate the result, but then would have to deal with
    ordering the annotations and make sure they don't overlap. We don't
    so make sure to keep these potential sources of error in mind.

    Note that the location indexing starts w/ 1, i.e. the index of the
    first position in the annotation is 1. Python starts indexing w/ 0.
    Don't worry, care has been taken. Note further that if the sequence is
    shorter than the annotation, a "truncated" subsequence is returned w/o
    any warning/ error messages.

    Example:

    from zoo.parse import location_todict, location_extract

    seq = 'ACTACACCCCGTGTAAA'
    loc = location_todict(['[1:10](-)', '[15:>30](?)'])  # len(seq) < 30
    assert(next(next(location_extract(loc, seq))) == 'AAA'
    '''
    for i in loc:
        yield seq[i['start']-1:i['end']]


