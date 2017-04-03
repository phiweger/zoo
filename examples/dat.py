'''shell
dat cp data/tests/cell_a.json tmp/dump.json
dat share send/
dat clone <link> receive/
'''


import json
import random


l = []
with open('tmp/dump.json', 'r+') as file:
    for line in file:
        l.append(
            json.loads(line.strip())
            )

# modify
for i in l:
    del i['dangerous']
    i['random'] = random.uniform(0, 1)


# dump again
with open('send/dump.json', 'w+') as file:
    for i in l:
        file.write(json.dumps(i))
        file.write('\n')


'''
dat pull receive/
# changes updated
'''
