rm -rf /tmp/mongorestore && mkdir /tmp/mongorestore
TMP_DIR=/tmp/mongorestore
# rm -rf $TMP_DIR && mkdir $TMP_DIR
if [[ $1 =~ \.tar$ ]];
then
	#FILENAME=$(echo $1 | sed 's/.*\///')
	FILENAME=$2"/"
	mkdir $TMP_DIR
	echo "Data will be extracted into :"$TMP_DIR
	tar -xvf $1 -C $TMP_DIR
    ls $TMP_DIR
else
	FILENAME=$(echo $1 | sed 's/.*\///')
	cp $1 $TMP_DIR$FILENAME
fi

docker run -it --rm --net net_mongodb \
-v $TMP_DIR:/tmp \
mongo bash -c \
'mongorestore --drop -v --host mongo:27017 --db '$2' /tmp/'$FILENAME
# rm -rf $TMP_DIR
rm -rf /tmp/mongorestore
