from copy import deepcopy
# from uuid import uuid4
from zoo.utils import deep_set, deep_get
from zoo.parse import location_tostr


def seqrecord2jsondict(seqrecord):
    '''In Bio.SeqRecord, out dict.

    TODO: jsondict2seqrecord()

    "jsondict" is a (nested) dictionary that is JSON serializable.

    Example:

    \b
    with open('KX882764.gb', 'r+') as file:
        a = SeqIO.read(file, 'genbank')
    print(json.dumps(seqrecord2json(a), indent=2))
    '''
    gb = deepcopy(seqrecord).__dict__

    # extract elements that are not JSON serializable:
    # features
    feat = gb.pop('features')
    gb['features'] = []
    for i in feat:
        i.location = next(location_tostr([i.location]))  # change in place
        gb['features'].append(i.__dict__)

    # references
    ref = []
    while True:
        try:
            r = deep_get(gb, 'annotations.references').pop()
            del r.__dict__['location']  # redundant, this info is in "features"
            ref.append(r.__dict__)
        except IndexError:  # pop from empty list
            break
    deep_set(gb, 'annotations.references', ref)

    # sequence
    gb['sequence'] = str(gb.pop('_seq'))
    return gb


# def keymap():
#     '''Map one dictionary onto another.

#     Use case: Map a genbank file to a zoo schema.
#     '''
#     pass


# def create_document(seqrecord, schema_base, schema_annotation):
#     # Bio.SeqRecord.SeqRecord
#     gb = deepcopy(seqrecord)
#     r = deepcopy(schema_base)  # r .. record

#     '''
#     Metadata.
#     '''

#     r['_id'] = str(uuid4())
#     # lets start enforcing the schema structure
#     deep_set(
#         r, 'metadata.alt_id',
#         {'gb': gb.name, 'gb_version': gb.id.split('.')[1]},
#         force=True, replace=True)

#     # db references
#     for i in gb.dbxrefs:
#         if i:
#             k, v = i.split(':')
#             deep_get(r, 'metadata.alt_id').update({k.lower(): v})

#     # description
#     deep_set(r, 'metadata.description', gb.description, force=True)

#     # gb.annotations
#     l = []
#     for i in gb.annotations['references']:
#         ref = deepcopy(i.__dict__)
#         del ref['location']  # redundant
#         l.append(ref)

#     del gb.annotations['references']
#     deep_set(r, 'metadata.references', l, force=True)
#     deep_set(r, 'metadata.dump', gb.annotations, force=True)
#     deep_set(
#         r, 'relative.taxonomy.dump',
#         r['metadata']['dump'].pop('taxonomy'), force=True)

#     '''
#     Annotations and a little metadata (feature of type "source").
#     '''

#     for f in gb.features:
#         if f.type != 'source':
#             a = deepcopy(schema_annotation)  # a .. annotation
#             a['source'] = 'gb'
#             a['type'] = f.type
#             a['location'] = location_tostr(f.location)

#             for k, v in {
#                 'name': 'product',
#                 'id': 'protein_id'
#             }.items():
#                 try:
#                     a[k] = f.qualifiers[v][0]
#                 except KeyError:
#                     pass
#             deep_get(r, 'derivative.annotation').append(a)
#         else:
#             for k, v in f.qualifiers.items():
#                 deep_set(r, 'metadata.source.' + k, v[0], force=True)
#     # db.zika.insert_one(r)
#     # db.zika.find_one()

#     '''
#     Sequence.
#     '''

#     deep_set(r, 'sequence', str(gb.seq))

#     return(r)


# for i in seqrecord:
#     create_document(i, base, annotation)