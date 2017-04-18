// node jsonpatch.js --db test --cell a --file diff_cell_a.json
// for line in nd-diff
//     get document by id
//     apply patch
//     replace document (or selectively modify keys)


var argv = require("minimist")(process.argv.slice(2));
var jsondiffpatch = require("jsondiffpatch");
var fs = require("fs");
var MongoClient = require("mongodb").MongoClient;
var mongoUtil = require( "mongoutils" );


// function get_collection(err, db) {
//     if (err) {
//         throw err;
//     }
//     // c = db.collection(argv.cell);
//     c = db.collection("a");
//     _id = "3";
//     c.find({"_id": _id}).limit(1).next(handle_document);
//     db.close()
// }


function handleDocument(err, doc){
    console.log(JSON.stringify(doc));
}

function closeConnection() {
    var db = mongoUtil.getDb();
    db.close()
}


// var c;
// var db = MongoClient.connect(
//     // "mongodb://" + argv.client + "/" + argv.db, 
//     "mongodb://" + "localhost:27017" + "/" + "test", 
//     // localhost == 127.0.0.1
//     // "mongodb://localhost:27017/test"
//     get_collection
//     );

// var fs = require("fs");
// var readStream = fs.createReadStream("diff.json");
// readStream.pipe(process.stdout);
// readStream.forEach(process.stdout);





function pump() {
    var pos;
    while ((pos = buf.indexOf("\n")) >= 0) { // keep going while there"s a newline somewhere in the buffer
        if (pos == 0) { // if there"s more than one newline in a row, the buffer will now start with a newline
            buf = buf.slice(1); // discard it
            continue; // so that the next iteration will start with data
        }
        processLine(buf.slice(0,pos)); // hand off the line
        buf = buf.slice(pos+1); // and slice the processed data off the buffer
    }
}


function processLine(line) { // here"s where we do something with a line
    if (line.length > 0) { // ignore empty lines
        var delta = JSON.parse(line); // parse the JSON
        var _id = delta._id
        delete delta._id
        console.log(delta, _id); // do something with the data here!

        var db = mongoUtil.getDb();
        db.collection("a").find({"_id": _id}).limit(1).next(handleDocument);
    }
}



// stream JSON file
// stackoverflow, 11874096
var stream = fs.createReadStream("diff.json", {flags: "r", encoding: "utf-8"});
var buf = "";


// reuse MongoDB connection
// stackoverflow, 24621940
mongoUtil.connectToServer( "test", function( err ) {

    stream.on("data", function(d) {
        buf += d.toString(); // when data is read, stash it in a string buffer
        pump(); // then process the buffer
        closeConnection();
    });
});









