// const express = require('express');
// const multer = require('multer');
// const mongoose = require('mongoose');
// const bodyParser = require('body-parser');

// const app = express();
// const port = process.env.PORT || 3000; // Set port for deployment

// // Configure Multer for file uploads
// const upload = multer({ dest: 'uploads/', fileFilter: (req, file, cb) => {
//   if (!file.originalname.match(/\.(fastq\.gz)$/)) {
//     return cb(new Error('Only FASTQ.gz files are allowed!'), false);
//   }
//   cb(null, true);
// } });
//  // Change 'uploads/' if needed

// // Connect to MongoDB
// mongoose.connect(`mongodb+srv://aryanb349:G00GLE34.Aryan@cluster0.qpjxkfz.mongodb.net/?retryWrites=true&w=majority&appName=Cluster0`, {
//   useNewUrlParser: true,
//   useUnifiedTopology: true
// })
//   .then(() => console.log('Connected to MongoDB'))
//   .catch(err => console.error('Error connecting to MongoDB:', err));

// // Define Mongoose schema for genome data
// const genomeSchema = new mongoose.Schema({
//   name: String,
//   fastqData: String, // Store the entire FASTA data as a string for simplicity
//   fastqcReport: String // Store the FASTQC report as a string
// });

// const uploadFile = require('./')

// const Genome = mongoose.model('Genome', genomeSchema);

// // Configure body-parser for form data access
// app.use(bodyParser.urlencoded({ extended: false }));

// // Route for uploading genome file
// aapp.post('/upload-genome', upload.single('genomeFile'), async (req, res) => {
//   try {
//     const { filename, originalname } = req.file;

//     // Read the uploaded file as a buffer
//     const fileBuffer = await fs.promises.readFile(`uploads/${filename}`);

//     // Decompress the buffer using 'zlib'
//     const zlib = require('zlib');
//     const decompressedData = zlib.gunzipSync(fileBuffer).toString('utf-8');

//     const genome = new Genome({
//       name: originalname,
//       fastaData: decompressedData, // Not suitable for FASTQ
//       fastqcReport: '' // Initialize empty report field
//     });

//     await genome.save();

//     // Call the Bash script to run FASTQC (see below)
//     const fastqcReport = await execSync(`bash fastqc_pipeline.sh ${filename}`); // Modify command if needed

//     // Update the genome data with the report
//     await Genome.updateOne({ _id: genome._id }, { fastqcReport: fastqcReport.toString() });

//     res.send('Genome uploaded and processed successfully!');
//   } catch (err) {
//     console.error('Error uploading or processing genome:', err);
//     res.status(500).send('Error processing genome. Please try again.');
//   }
// });

// // Route to serve FASTQC report (assuming the report is stored in the database)
// app.get('/report/:genomeId', async (req, res) => {
//   try {
//     const genome = await Genome.findById(req.params.genomeId);
//     if (!genome) {
//       return res.status(404).send('Genome not found.');
//     }

//     res.send(genome.fastqcReport); // Send the report data
//   } catch (err) {
//     console.error('Error fetching FASTQC report:', err);
//     res.status(500).send('Error retrieving report. Please try again.');
//   }
// });

// app.listen(port, () => console.log(`Server listening on port ${port}`));

const express = require('express');
const multer = require('multer');
const { exec } = require('child_process');
const fs = require('fs');
const path = require('path');

const app = express();
const port = 3000;

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

// Serve the MultiQC report
app.get('/multiqc_report', (req, res) => {
    // Send the MultiQC report as response
    const reportPath = path.join(__dirname, 'multiqc_reports', 'multiqc_report.html');
    res.sendFile(reportPath);
});

app.listen(port, () => {
    console.log(`Server is running at http://localhost:${port}`);
});


