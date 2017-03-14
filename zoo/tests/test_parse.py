from zoo.parse import parse_date


def test_parse_date():
    date = '1976/1/09'
    assert parse_date(date)['y'] == 1976
