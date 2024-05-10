const express = require('express');
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




// Multer configuration for file upload
const storage = multer.diskStorage({
    destination: function (req, file, cb) {
        cb(null, 'uploads/')
    },
    filename: function (req, file, cb) {
        cb(null, file.originalname)
    }
});
const upload = multer({ storage: storage });


// mongoose.connect(`mongodb://localhost:27017/omics_db`, { useNewUrlParser: true, useUnifiedTopology: true });


// Serve static files
app.use(express.static('public'));

app.post('/upload', upload.array('genomeFiles', 2), (req, res) => {
    // Paths to uploaded files
    const filePaths = req.files.map(file => file.path);

    // Run the Bash script
    exec(`bash fastqc.sh "${filePaths[0]}" "${filePaths[1]}"`, (error, stdout, stderr) => {
        if (error) {
            console.error(`exec error: ${error}`);
            res.status(500).send('Error occurred during processing');
            return;
        }
        // Read the generated report
        const reportPath = path.join(__dirname, 'multiqc_reports', 'multiqc_report.html');
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

// app.post('/upload', upload.array('genomeFiles', 2), (req, res) => {
//   // Paths to uploaded files
//   const filePaths = req.files.map(file => file.path);

//   // Extract base names and group files
//   const fileGroups = {};
//   req.files.forEach(file => {
//       const baseName = file.originalname.split('_')[0];
//       if (!fileGroups[baseName]) {
//           fileGroups[baseName] = [];
//       }
//       fileGroups[baseName].push(file.originalname);
//   });

//   // Convert to desired JSON format
//   const filePairs = Object.keys(fileGroups).map(baseName => {
//       const files = fileGroups[baseName];
//       return {
//           file_name: baseName,
//           file_1: files[0],
//           file_2: files[1]
//       };
//   });

//   // Get the most recent pair
//   const mostRecentPair = filePairs[filePairs.length - 1];

//   // Run the Bash script on the most recent pair
//   exec(`bash fastqc.sh "${mostRecentPair.file_1}" "${mostRecentPair.file_2}"`, (error, stdout, stderr) => {
//       if (error) {
//           console.error(`exec error: ${error}`);
//           res.status(500).send('Error occurred during processing');
//           return;
//       }
//       // Read the generated report
//       const reportPath = path.join(__dirname, 'multiqc_reports', 'multiqc_report.html');
//       fs.readFile(reportPath, 'utf8', (err, data) => {
//           if (err) {
//               console.error(`Error reading report: ${err}`);
//               res.status(500).send('Error occurred during processing');
//               return;
//           }
//           // Send the report as response
//           res.send(data);
//       });
//   });
// });

// Serve the MultiQC report
app.get('/multiqc_report', (req, res) => {
    // Send the MultiQC report as response
    const reportPath = path.join(__dirname, 'multiqc_reports', 'multiqc_report.html');
    res.sendFile(reportPath);
});

app.listen(port, () => {
    console.log(`Server is running at http://localhost:${port}`);
});




