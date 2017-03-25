import json


with open('survey.json', 'w+') as outjson:
    outjson.write(dumps(db.survey.find(), indent=4))
# fasta is 350 kb, JSON 359


with open('survey.json', 'w+') as output:
    for i in db.survey.find():
        output.write(json.dumps(i) + '\n')



'''
mongoexport --db zika --collection survey --out survey.json

{"_id":"86853586-5e9...
{"_id":"689e59b8-514...
{"_id":"6d9bff35-aab...

we can then do
with open...
    for line in input.readline()  # i.e. streaming
'''