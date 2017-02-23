rm -rf /tmp/mongodump && mkdir /tmp/mongodump

# docker run -it --rm --link mongo:mongo \
docker run -it --rm --net net_mongodb \
-v /tmp/mongodump:/tmp \
mongo bash -c \
'mongodump -v --host mongo:27017 --db '$1' --out=/tmp'
# docker run -it --rm --net net_mongodb mongo mongo --host mongo:27017

ls -ahl /tmp/mongodump
tar -cvf $2 -C /tmp/mongodump .
rm -rf /tmp/mongodump

