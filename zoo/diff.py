from deepdiff import DeepDiff
import json


'''
wrap this go lib?
https://github.com/benjamine/jsondiffpatch

native python
http://json-delta.readthedocs.io/en/latest/philosophy.html
https://github.com/ZoomerAnalytics/jsondiff

# go get github.com/yudai/gojsondiff/jp
jd a.json b.json | jq -c . | tr -d '^J' > .zoo/diff/someid.delta
# https://github.com/stedolan/jq/issues/215
jp .zoo/diff/someid.delta a.json
'''


def jd(collection, out, file):
    '''Diff a.json and b.json into nd-json.delta

    This function simply diffs 2 JSON files into a newline-delimited JSON
    file. This allows to concatenate diffs to a delta file, and a rollback
    will only need to apply all the changes from the bottom to the top of the
    delta file.

    There is a separate patch function ("jp"). The idea is to at some time
    in the future (written 2017-04-13) replace the external dependency for
    a native diff implementation.

    # go get github.com/yudai/gojsondiff
    jd -f delta a.json b.json | jq -c . | tr -d '^J' > .zoo/diff/someid.delta
    # https://github.com/stedolan/jq/issues/215
    jp .zoo/diff/someid.delta a.json

    or
    https://github.com/benjamine/jsondiffpatch
    and write some js code - load json, diff, write.
    or
    stackoverflow, 28810667
    https://www.npmjs.com/package/deep-diff
    '''
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


def jp():
    '''Apply a patch (.delta) to a JSON file.'''
    pass
