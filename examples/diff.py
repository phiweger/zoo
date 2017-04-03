from deepdiff import DeepDiff
import json
from pprint import pprint
from pymongo import MongoClient
from zoo import get_data


fp_a = get_data('tests/cell_a.json')
fp_b = get_data('tests/cell_b.json')
fp_diff = 'diff.json'


c = MongoClient('localhost:27017')['diff']['mock']
with open(fp_a, 'r+') as file:
    for line in file:
        c.insert_one(
            json.loads(line.strip())
            )


with open(fp_b, 'r+') as file, open(fp_diff, 'w+') as out:
    for line in file:
        old = json.loads(line.strip())
        _id = old['_id']
        new = c.find_one({'_id': _id})
        diff = DeepDiff(
            new, old,  # note order of new, old
            ignore_order=True,
            verbose_level=2)
        if diff:
            diff['_id'] = _id
            out.write(diff.json)
            out.write('\n')
        pprint(diff, indent=2)

