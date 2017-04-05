from zoo.align import encode_gaps, decode_gaps, hash_seq


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


def test_encode_decode():
    seq = 'ACTGAGGTAGGTAA'
    seq_gap = 'A--CTGA---GGTAGGT-AA'
    assert decode_gaps(seq, encode_gaps(seq_gap)) == seq_gap


def test_hash_seq(hash='md5'):
    assert hash_seq(['AAAA', 'ACTG']) == '3b0d5a673fb29e5201e4587a35e2576d'
    assert hash_seq(['ACTG', 'AAAA']) == hash_seq(['AAAA', 'ACTG'])

