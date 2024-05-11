const express = require('express');
const bodyParser = require('body-parser');
const crypto = require("crypto");
const multer = require('multer');
const { exec } = require('child_process');
const fs = require('fs');
const path = require('path');
const mongoose = require("mongoose");
const app = express();
const port = 3000;

// add env vars from env for mongodb
const mongo_user = process.env.mongo_user;
const mongo_pass = process.env.mongo_pass;

// Add middleware for parsing application/x-www-form-urlencoded
app.use(bodyParser.urlencoded({ extended: true }));

// Add middleware for parsing application/json
app.use(bodyParser.json());

// Multer configuration for file upload
const storage = multer.diskStorage({
    destination: function (req, file, cb) {
        const token = req.body.token ; // Extract token from request body
        if (!token) {
            // Handle missing token error (e.g., return error or default path)
            return cb(new Error('Missing token in request body'));
        }
        const uploadDir = path.join('uploads', token); // Create upload directory with token
        if (!fs.existsSync(uploadDir)) {
            fs.mkdirSync(uploadDir, { recursive: true }); // Create directory recursively if not exists
        }
        cb(null, uploadDir);
    },
    filename: function (req, file, cb) {
        cb(null, file.originalname)
    }
});
const upload = multer({ storage: storage });

app.use(express.static('public'));

app.post('/upload', upload.array('genomeFiles', 2), (req, res) => {
    // Paths to uploaded files
    const token = req.body.token; // Extract token from request body
    const filePaths = req.files.map(file => file.path);

    // Run the Bash script
    exec(`bash fastqc.sh "${filePaths[0]}" "${filePaths[1]}"`, (error, stdout, stderr) => {
        if (error) {
            console.error(`exec error: ${error}`);
            res.status(500).send('Error occurred during processing');
            return;
        }
        const reportDir = path.join('multiqc_reports', token); // Create report directory with token
        if (!fs.existsSync(reportDir)) {
            fs.mkdirSync(reportDir, { recursive: true }); // Create directory recursively if not exists
        }

        // Read the generated report
        const reportPath = path.join(reportDir, 'multiqc_report.html');
        fs.readFile(reportPath, 'utf8', (err, data) => {
            if (err) {
                console.error(`Error reading report: ${err}`);
                res.status(500).send('Error occurred during processing');
                return;
            }
            // Send the report as response
            res.send(data);
        });
    });
});

// Serve the MultiQC report
app.get('/multiqc_report/:token', (req, res) => {
    const token = req.params.token;
    const reportDir = path.join('multiqc_reports', token);
    const reportPath = path.join(reportDir, 'multiqc_report.html');
    res.sendFile(reportPath);
});

app.listen(port, () => {
    console.log(`Server is running at http://localhost:${port}`);
});
