const express = require('express');
const crypto = require("crypto");
const multer = require('multer');
const { exec } = require('child_process');
const fs = require('fs');
const path = require('path');
const app = express();
const port = 3000;

// Define storage for uploaded files
const storage = multer.diskStorage({
    destination: function (req, file, cb) {
        cb(null, 'uploads/')
    },
    filename: function (req, file, cb) {
        cb(null, file.originalname)
    }
});
const upload = multer({ storage: storage });
console.log("uploads")
// Serve static files
app.use(express.static('public'));

// Endpoint for file upload
app.post('/upload', upload.array('genomeFiles', 2), (req, res) => {
    // Get the token from the request
    const token = req.body.token;
    // Paths to uploaded files
    const filePaths = req.files.map(file => file.path);

    // Move uploaded files to a directory with the token as name
    const uploadDir = path.join(__dirname, 'uploads', token);
    fs.mkdirSync(uploadDir, { recursive: true });
    for (let i = 0; i < filePaths.length; i++) {
        const destPath = path.join(uploadDir, path.basename(filePaths[i]));
        fs.renameSync(filePaths[i], destPath);
    }
    console.log("uploadDir")
    // Run the Bash script
    exec(`bash fastqc.sh "${uploadDir}"`, (error, stdout, stderr) => {
        if (error) {
            console.error(`exec error: ${error}`);
            res.status(500).send('Error occurred during processing');
            return;
        }
    console.log("exec")
        // Read the generated report
        const reportPath = path.join(__dirname, 'multiqc_reports', token, 'multiqc_report.html');
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

// Endpoint for serving MultiQC report
app.get('/multiqc_report/:token', (req, res) => {
    const token = req.params.token;
    // Send the MultiQC report as response
    const reportPath = path.join(__dirname, 'multiqc_reports', token, 'multiqc_report.html');
    res.sendFile(reportPath);
});

app.listen(port, () => {
    console.log(`Server is running at http://localhost:${port}`);
});
