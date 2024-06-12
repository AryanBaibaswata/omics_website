const express = require('express');
const multer = require('multer');
const { exec } = require('child_process');
const path = require('path');
const fs = require('fs');
const app = express();
require('dotenv').config(); 
require('./config/passport')(passport);
const PORT = process.env.PORT || 3000;


app.use(express.static(path.join(__dirname, 'public')));
