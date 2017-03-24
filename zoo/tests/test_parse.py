from zoo.parse import parse_date
from zoo.parse import location_todict, location_tostr
from Bio.SeqFeature import FeatureLocation, AfterPosition


def test_parse_date():
    date = '1976/1/09'
    assert parse_date(date)['y'] == 1976


def test_location_tostr():
    f1 = FeatureLocation(5, 10, strand=-1)
    f2 = FeatureLocation(20, AfterPosition(30), strand=0)
    combined = f1 + f2
    assert location_tostr(combined) == ['[5:10](-)', '[20:>30](?)']


def test_location_todict():
    loc = ['[5:10](-)', '[20:>30](?)']
    assert location_todict(loc)[1]['partial'] is True
