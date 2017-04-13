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


def json_diff(collection, out, file):
    '''Diff document from MongoDB with line in nd-JSON, out diff.delta.

    This function simply diffs documents into a newline-delimited JSON
    file. This allows to concatenate diffs to a delta file, and a rollback
    will only need to apply all the changes from the bottom to the top of the
    delta file.

    There is a separate patch function ("json_patch"). The idea is to at some
    time in the future (written 2017-04-13) replace the external dependency for
    a native diff implementation.

    # go get github.com/yudai/gojsondiff
    jd -f delta a.json b.json | jq -c . | tr -d '^J' > .zoo/diff/someid.delta
    # https://github.com/stedolan/jq/issues/215
    jp .zoo/diff/someid.delta a.json

    # https://github.com/stefankoegl/python-json-patch

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

    # from subprocess import call
    # Note that you have to specify path to script
    # success = call(["node", "/Users/pi/repos/zoo/wip/hello.js"])
    # 1. get line new
    # 2. get diff to old in database
    # 3. export diff as json

    # then checkout
    # iterate over diffs and path them till id of diff reached
    # git log shows ids of commits (incl. diffs) w/ meassage and datetime


def json_patch():
    '''Apply a patch (.delta) to a JSON file.'''
    pass
