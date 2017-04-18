// node jsonpatch.js --db test --cell a --file diff_cell_a.json
// for line in nd-diff
//     get document by id
//     apply patch
//     replace document (or selectively modify keys)



// node ~/repos/zoo/zoo/scripts/jsonpatch.js --client localhost:27017 --db test --cell a --file diff.json


var argv = require("minimist")(process.argv.slice(2));
var jsondiffpatch = require("jsondiffpatch");
var fs = require("fs");
var MongoClient = require("mongodb").MongoClient;
var mongoutils = require( "mongoutils" );
// var sync = require( "async" );


function handleDocument(err, doc){
    // console.log(JSON.stringify(doc));
    // console.log(this._id);
    // console.log(this.delta);
    try{
        var patch = jsondiffpatch.patch(doc, this.delta);
    }
    catch (err){
        console.log('The delta in inconsistent with the document state.')
        process.exit()
    }
    // console.log(patch);
    var db = mongoutils.getDb();
    db.collection(argv.cell)
        .update({"_id": this._id}, {patch}, function (err, result){
            counter -= 1;
            if (counter == 0){
                // stackoverflow, 8373905
                // pump() increases the counter for each line of JSON
                // every update decreases the counter
                closeConnection();
            }
            // console.log("counter:", counter);
        });
    //return cb()
    // insertDocument(this._id, patch)
    // callback()
}


function closeConnection() {
    var db = mongoutils.getDb();
    db.close()
}


function pump() {
    var pos;
    while ((pos = buf.indexOf("\n")) >= 0) { 
    // keep going while there"s a newline somewhere in the buffer
        if (pos == 0) { // if there"s more than one newline in a row, the buffer will now start with a newline
            buf = buf.slice(1); // discard it
            continue; // so that the next iteration will start with data
        }
        processLine(buf.slice(0, pos)); // hand off the line
        buf = buf.slice(pos + 1); // and slice the processed data off the buffer
        counter += 1;
    }
}


// here's where we do something with a line
function processLine(line) {
    if (line.length > 0) { // ignore empty lines
        
        var delta = JSON.parse(line); // parse the JSON
        var _id = delta._id
        delete delta._id
        // console.log(delta, _id); // do something with the data here!

        var db = mongoutils.getDb();
        db.collection(argv.cell)
            .find({"_id": _id})
            .limit(1)
            .next(handleDocument.bind({
                delta: delta,
                _id: _id
            }));
            // bind() to pass variables to callback functions
            // http://samwize.com/2013/09/01/how-you-can-pass-a-variable-\
            // into-callback-function-in-node-dot-js/
    }
}


// stream JSON file
// stackoverflow, 11874096
// there's a dedicated module: https://github.com/dominictarr/JSONStream
var stream = fs.createReadStream(argv.file, {flags: "r", encoding: "utf-8"});
var buf = "";
var counter = 0


// reuse MongoDB connection
// stackoverflow, 24621940
mongoutils.connectToServer(argv.client, argv.db, function( err ) {

    stream.on("data", function(d) {
        buf += d.toString(); // when data is read, stash it in a string buffer
        pump(); // then process the buffer
        // closeConnection();
    });
});


