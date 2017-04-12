from deepdiff import DeepDiff
import json


def diffutil(collection, out, file):
    with open(out, 'w+') as outfile:
        for line in file:
            old = json.loads(line.strip())
            _id = old['_id']
            new = collection.find_one({'_id': _id})
            diff = DeepDiff(
                new, old,  # note order of new, old
                ignore_order=True,
                verbose_level=2)
            if diff:
                diff['_id'] = _id
                outfile.write(diff.json)
                outfile.write('\n')
