const mongoose = require('mongoose');
const gridfsStream = require('gridfs-stream');
const { GridFSBucket } = require('mongodb');

// Connect to MongoDB
const conn = mongoose.createConnection('mongodb://localhost:27017/yourDB', {
    useNewUrlParser: true,
    useUnifiedTopology: true
});

let gfs, gridFSBucket;

conn.once('open', () => {
    gridFSBucket = new GridFSBucket(conn.db, {
        bucketName: 'uploads'
    });
    gfs = gridfsStream(conn.db, mongoose.mongo);
    gfs.collection('uploads');
});

module.exports = { gfs, gridFSBucket };
