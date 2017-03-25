from zoo.parse import parse_date
from zoo.parse import location_todict, location_tostr, location_extract
from Bio.SeqFeature import FeatureLocation, AfterPosition
import pytest


def test_parse_date():
    date = '1976/1/09'
    assert parse_date(date)['y'] == 1976


def test_location_tostr_from_featloc():
    '''Input is an iterable of Biopython FeatureLocation (can be Compound).'''
    f1 = FeatureLocation(5, 10, strand=-1)
    f2 = FeatureLocation(20, AfterPosition(30), strand=0)
    combined = f1 + f2
    loc = [f1, f2, combined]
    gen = location_tostr(loc)
    assert next(gen) == '[5:10](-)'
    assert next(gen) == '[20:>30](?)'
    assert next(gen) == ['[5:10](-)', '[20:>30](?)']
    with pytest.raises(StopIteration):
        next(gen)


def test_location_tostr_fromdict():
    ''''Input is an iterable of dicts such as from location_todict().'''
    loc = ['[5:10](-)', '[20:>30](?)']
    gen = location_todict(loc)
    assert next(location_tostr(gen)) == loc[0]


def test_location_todict():
    loc = ['[5:10](-)', '[20:>30](?)']
    gen = location_todict(loc)
    assert next(gen)['start'] == 5
    assert {v for k, v in next(gen).items()} == {20, 30, 1, 0, '?'}


def test_location_extract():
    seq = 'ACTACACCCCGTGTAAA'
    loc = location_todict(['[1:10](-)', '[15:>30](?)'])  # len(seq) < 30
    gen = location_extract(loc, seq)
    assert next(gen) == 'ACTACACCCC'
    assert next(gen) == 'AAA'  # "trunctated" because annptation extends beyond
