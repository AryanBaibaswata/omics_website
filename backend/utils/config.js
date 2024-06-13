require('dotenv').config;

const port = process.env.PORT || 3000;
const mongodb_URI = process.env.mongodb_URI;


module.exports = {port, mongodb_URI}