from zoo.align import encode_gaps, decode_gaps


def test_encode_gaps():
    seq = 'A--CTGA---GGTAGGT-AA'
    assert encode_gaps(seq) == {1: 2, 7: 3, 17: 1}


def test_decode_gaps():
    '''
    Appropriately handly upper and lower case, because e.g. Mafft alignments
    are returned (by default) in lower case letters.
    '''
    seq = 'AcTGAggtAGGTAA'
    gapdict = {1: 2, 7: 3, 17: 1}
    result = decode_gaps(seq, gapdict, uppercase=False)
    assert result == 'A--cTGA---ggtAGGT-AA'
