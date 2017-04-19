var MongoClient = require('mongodb').MongoClient;
var _db;

module.exports = {

  connectToServer: function (client, cell, callback){
    MongoClient.connect("mongodb://" + client + "/" + cell, function(err, db) {
      _db = db;
      return callback(err);
    } );
  },

  getDb: function(){
    return _db;
  }
};