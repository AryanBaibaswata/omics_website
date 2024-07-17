const express = require('express');
const multer = require('multer');
const { exec } = require('child_process');
const path = require('path');
const fs = require('fs');
const app = express();
const PORT = process.env.PORT || 3000;

const UPLOAD_FOLDER = path.join(__dirname, 'uploads');

// Set up multer for file uploads
const storage = multer.diskStorage({
    destination: (req, file, cb) => {
        const date = new Date().toISOString().replace(/:/g, '-');
        req.uploadDate = date; // Store the date for this upload session
        const uploadPath = path.join(UPLOAD_FOLDER, date, 'files');
        fs.mkdirSync(uploadPath, { recursive: true });
        fs.mkdirSync(path.join(UPLOAD_FOLDER, date, 'fastqc_output'), { recursive: true });
        cb(null, uploadPath);
    },
    filename: (req, file, cb) => {
        cb(null, file.originalname);
    }
});

const upload = multer({ storage });

// Middleware to serve static files
app.use(express.static(path.join(__dirname, 'public')));

// Helper functions
function validateFiles(files, refGenome) {
    if (files.length % 2 !== 0) {
        throw new Error('Please upload files in pairs.');
    }
    if (!refGenome || (!refGenome.originalname.endsWith('.fasta') && !refGenome.originalname.endsWith('.fa'))) {
        throw new Error('Please upload a .fasta or .fa file as reference genome.');
    }
}

function processFiles(files, basedir) {
    const sampleFiles = [];
    for (let i = 0; i < files.length; i += 2) {
        const file1 = files[i].filename;
        const file2 = files[i + 1].filename;
        
        if (file1.replace('_1.fastq.gz', '') === file2.replace('_2.fastq.gz', '')) {
            sampleFiles.push(file1.replace('_1.fastq.gz', ''));
        } else if (file1.replace('_R1.fastq.gz', '') === file2.replace('_R2.fastq.gz', '')) {
            const sampleName = file1.replace('_R1.fastq.gz', '');
            sampleFiles.push(sampleName);
            renameFiles(file1, file2, sampleName, basedir);
        } else {
            throw new Error('File pairs do not match.');
        }
    }
    return sampleFiles;
}

function renameFiles(file1, file2, sampleName, basedir) {
    const oldPath1 = path.join(basedir, file1);
    const oldPath2 = path.join(basedir, file2);
    const newPath1 = path.join(basedir, `${sampleName}_1.fastq.gz`);
    const newPath2 = path.join(basedir, `${sampleName}_2.fastq.gz`);
    fs.renameSync(oldPath1, newPath1);
    fs.renameSync(oldPath2, newPath2);
}

function generatePipelineScript(basedir, samplesList, GENOMEIDX1, GENOMEIDX) {
    // ... (keep the existing pipeline script generation logic)
}

// Route to handle file uploads and trigger the pipeline
app.post('/upload', upload.fields([
    { name: 'files', maxCount: 100 },
    { name: 'refgenome', maxCount: 1 }
]), (req, res) => {
    try {
        const files = req.files['files'].sort();
        const refGenome = req.files['refgenome'][0];
        validateFiles(files, refGenome);

        const basedir = path.join(UPLOAD_FOLDER, req.uploadDate);
        
        const genomeFilePath = path.join(basedir, 'files', refGenome.filename);
        const indexBasePath = path.join(basedir, path.basename(refGenome.filename, path.extname(refGenome.filename)));
        
        const bowtieBuildCmd = `bowtie2-build ${genomeFilePath} ${indexBasePath}`;

        exec(bowtieBuildCmd, (error, stdout, stderr) => {
            if (error) {
                console.error(`Error building Bowtie index: ${error}`);
                return res.status(500).send('Error building Bowtie index.');
            }

            const sampleFiles = processFiles(files, path.join(basedir, 'files'));
            const samplesList = sampleFiles.map(sample => `    "${sample}"`).join(' \\\n');
            
            const pipelineScriptContent = generatePipelineScript(basedir, samplesList, indexBasePath, genomeFilePath);

            const scriptPath = path.join(basedir, 'pipeline.sh');
            fs.writeFileSync(scriptPath, pipelineScriptContent);
            fs.chmodSync(scriptPath, '755');

            const progressFile = path.join(basedir, 'progress.txt');
            fs.writeFileSync(progressFile, ''); // Clear previous progress

            const child = exec(`bash ${scriptPath}`, { shell: '/bin/bash' });

            child.stdout.on('data', (data) => {
                console.log(`stdout: ${data}`);
                fs.appendFileSync(progressFile, data);
            });

            child.stderr.on('data', (data) => {
                console.error(`stderr: ${data}`);
                fs.appendFileSync(progressFile, data);
            });

            child.on('close', (code) => {
                console.log(`child process exited with code ${code}`);
                fs.appendFileSync(progressFile, `Process completed with code ${code}\n`);
            });

            res.json({ message: 'Pipeline execution started.' });
        });
    } catch (error) {
        res.status(400).json({ error: error.message });
    }
});

// ... (keep the existing routes for progress, download-progress, and file listing)

// Start the server
app.listen(PORT, () => {
    console.log(`Server is running on http://localhost:${PORT}`);
});