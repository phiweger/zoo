// node hello.js --db test --cell a --file cell_a.json

// works
// also check out:
// https://github.com/mafintosh/mongojs
// http://mongoosejs.com/
// http://mongoosejs.com/docs/index.html


// http://doduck.com/node-js-mongodb-hello-world-example/
var argv = require('minimist')(process.argv.slice(2));
var jsondiffpatch = require('jsondiffpatch');
var fs = require('fs');
var MongoClient = require('mongodb').MongoClient;


console.dir(argv);
// var name_db = 'test' // argv.db  // test
// var name_collection = 'a'//argv.cell  // a
// console.log(argv.db)
// console.log(argv.cell)


var c;
var db = MongoClient.connect(
    'mongodb://127.0.0.1:27017/' + argv.db, 
    function(err, db) {
        if(err)
            throw err;
        // console.log('connected to the mongoDB !');
        c = db.collection(argv.cell);
        // console.log(c.namespace)  // test.a
        // db.close();
        fs.readFileSync(argv.file).toString().split('\n').forEach(
            function (line) { // err, line
            //    console.error(err); // throw err?
            // process.exit()
            var l = JSON.parse(line);
            // console.log(l._id);

            c.findOne({'_id': l._id}, function(err, doc) {
                if(err) {
                    throw err;
                    }
                var delta = jsondiffpatch.diff(doc, l);
                delta._id = l._id
                // console.log(delta)
                result = JSON.stringify(delta)
                fs.appendFileSync('./commit_id.json', result + '\n');
            });
        });
        db.close()
    });

