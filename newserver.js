const path = require('path');
const express = require('express');
const bodyParser = require('body-parser');
const multer = require('multer');
const { exec } = require('child_process');
const fs = require('fs');
const app = express();
const port = 3000;
const date = new Date().toISOString().replace(/:/g, '-');
const uploadPath = path.join(__dirname, `uploads/${date}/files`);
// Create a storage object with the desired destination
const storage = multer.diskStorage({
    destination: (req, file, cb) => {

        fs.mkdirSync(uploadPath, { recursive: true });
        cb(null, uploadPath);
    },
    filename: (req, file, cb) => {
        cb(null, file.originalname);
    }
});

// Create the multer middleware using the storage object
const upload = multer({ storage: storage });

app.use(express.static('public'));

app.post('/upload', upload.array('genomeFiles', 20), (req, res) => {
    const filePaths = req.files.map(file => file.path);
    exec(`bash fastqc.sh ${filePaths.join(' ')} ${date}`, (error, stdout, stderr) => {
        if (error) {
            console.error(`exec error: ${error}`);
            res.status(500).send('Error occurred during processing');
            return;
        }
        const reportPath = path.join(__dirname, `multiqc_reports/${date}`, 'multiqc_report.html');
        fs.readFile(reportPath, 'utf8', (err, data) => {
            if (err) {
                console.error(`Error reading report: ${err}`);
                res.status(500).send('Error occurred during processing');
                return;
            }
            res.send(data);
        });
    });
});

app.get('/multiqc_report', (req, res) => {
    // Send the MultiQC report as response
    const reportPath = path.join(__dirname, `multiqc_reports/${date}`, 'multiqc_report.html');
    res.sendFile(reportPath);
});

app.listen(port, () => {
    console.log(`Server is running on port ${port}`);
});