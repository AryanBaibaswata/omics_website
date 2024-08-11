const express = require('express');
const multer = require('multer');
const { exec } = require('child_process');
const path = require('path');
const fs = require('fs');
const zlib = require('node:zlib')
const json = require('json');
const util = require('util');
const execPromise = util.promisify(exec);
const renamePromise = util.promisify(fs.rename);
// // const passport = require('passport');
// // const LocalStrategy = require('passport-local').Strategy;
// // const session = require('express-session');
// // const redis = require("connect-redis")(session);
// // const bcrypt = require('bcrypt');
const mongoose = require('mongoose');
const authRoutes = require('./controllers/auth');
const authMiddleware = require('./middleware/auth');
require('dotenv').config();
const app = express();
const PORT = process.env.PORT || 3000;
// const handlebars = require('express-handlebars');
// const handlebars = require('express-handlebars');
// Set up multer for file uploads

mongoose.connect(process.env.MONGODB_URI, {
    useNewUrlParser: true,
    useUnifiedTopology: true,
});

app.use(express.json());
app.use('/auth', authRoutes);



function generateUniqueFolderName() {
    const currentDate = new Date();
    const dateIST = new Date(currentDate.getTime());
    const options = {
        year: 'numeric',
        month: '2-digit',
        day: '2-digit',
        hour: '2-digit',
        minute: '2-digit',
        second: '2-digit',
        hour12: false,
        timeZone: 'Asia/Kolkata'
    };
    const formatter = new Intl.DateTimeFormat('en-GB', options);
    const formattedDate = formatter.format(dateIST);
    const formattedIST = formattedDate.replace(/\/|:|,/g, '-').replace(' ', 'T');
    return formattedIST;
}

app.use((req, res, next) => {
    req.uploadFolder = generateUniqueFolderName();
    next();
});


const storage = multer.diskStorage({
    destination: (req, file, cb) => {
        const uploadPath = `uploads/${req.uploadFolder}/files`;
        fs.mkdirSync(uploadPath, { recursive: true });
        fs.mkdirSync(`uploads/${req.uploadFolder}/fastqc_output`, { recursive: true });
        cb(null, uploadPath);
    },
    filename: (req, file, cb) => {
        cb(null, file.originalname);
    }
});


const upload = multer({ storage }).fields([
    { name: 'files' }, // Adjust maxCount as needed
    { name: 'refgenome', maxCount: 1 }
]);
// console.log("read part till storage")

// Middleware to serve static files
app.use(express.static(path.join(__dirname, 'public')));
// console.log("read till middleware")
// Route to handle file uploads and trigger the pipeline
app.post('/upload', upload, async (req, res) => {
    const files = req.files['files'];
    const genomeFiles = req.files['refgenome'];
    if (!genomeFiles || genomeFiles.length === 0) {
        return res.status(400).send('Please upload a .fasta or .fa file.');
    }
    try {
        // Add this new preprocessing step
        await Promise.all(files.map(async (file) => {
            let filePath = path.join(`uploads/${req.uploadFolder}/files`, file.filename);

            if (file.filename.endsWith('.fq')) {
                // Rename .fq to .fastq and gzip
                let newFilename = file.filename.slice(0, -3) + '.fastq';
                let newFilePath = path.join(`uploads/${req.uploadFolder}/files`, newFilename);
                await renamePromise(filePath, newFilePath);
                await execPromise(`gzip ${newFilePath}`);
                file.filename = newFilename + '.gz';
                console.log("renaming fq");
            } else if (file.filename.endsWith('.fastq')) {
                // Just gzip
                await execPromise(`gzip ${filePath}`);
                file.filename += '.gz';
                console.log('renaming .fastq');
            } else if (file.filename.endsWith('.fq.gz')) {
                // Gunzip, rename to .fastq, then gzip again
                console.log("renaming .fq.gz")
                await execPromise(`gunzip ${filePath}`);
                let unzippedPath = filePath.slice(0, -3);
                let newFilename = file.filename.slice(0, -6) + '.fastq';
                let newFilePath = path.join(`uploads/${req.uploadFolder}/files`, newFilename);
                await renamePromise(unzippedPath, newFilePath);
                await execPromise(`gzip ${newFilePath}`);
                file.filename = newFilename + '.gz';
            }
        }))
    } catch (err) {
        console.error('Error during file preprocessing:', err);
        return res.status(500).send('Error during file preprocessing.');
    }



    const genomeFile = genomeFiles[0];
    let genomeFilePath = path.join(`uploads/${req.uploadFolder}/files`, genomeFile.filename);
    //testing for fq and fq.gz files 


    // Replace spaces with underscores in genomeFile.filename
    if (genomeFile.filename.includes(' ')) {
        const newFilename = genomeFile.filename.replace(/\s+/g, "_");
        const newFilePath = path.join(`uploads/${req.uploadFolder}/files`, newFilename);
        fs.renameSync(genomeFilePath, newFilePath);
        genomeFilePath = newFilePath;
        genomeFile.filename = newFilename;
    }

    if (genomeFile.filename.endsWith('.fas')) {
        const newFilename = genomeFile.filename.slice(0, -4) + '.fasta';
        const newFilePath = path.join(`uploads/${req.uploadFolder}/files`, newFilename);
        fs.renameSync(genomeFilePath, newFilePath);
        genomeFilePath = newFilePath;
        genomeFile.filename = newFilename;
    }

    console.log(genomeFilePath);
    if (files.length % 2 !== 0) {
        return res.status(400).send('Please upload files in pairs.');
    }

    const genomeFilename = path.basename(genomeFile.filename, path.extname(genomeFile.filename));
    console.log(genomeFilename);

    // Replace spaces with underscores in the directory name
    const genomeDirName = genomeFilename.replace(/\s+/g, "_");
    fs.mkdirSync(`uploads/${req.uploadFolder}/${genomeDirName}`);

    // Replace spaces with underscores in indexBasePath
    const indexBasePath = path.join(`uploads/${req.uploadFolder}/${genomeDirName}`, genomeDirName);
    console.log("indexbasepath:", indexBasePath);

    const bowtieBuildCmd = `bowtie2-build ${genomeFilePath} ${indexBasePath}`;
    exec(bowtieBuildCmd, (error, stdout, stderr) => {
        console.log("began building")
        if (error) {
            console.error(`Error building Bowtie index: ${error}`);
            return res.status(500).send('Error building Bowtie index.');
        }

        console.log("reading files")
        const sampleFiles = [];
        console.log(sampleFiles)
        for (let i = 0; i < files.length; i += 2) {
            const file1 = files[i].filename;
            const file2 = files[i + 1].filename;
            console.log(file1, file2)
            if (file1.replace('_1.fastq.gz', '') === file2.replace('_2.fastq.gz', '')) {
                sampleFiles.push(file1.replace('_1.fastq.gz', ''));
            } else if (file2.replace('_1.fastq.gz', '') === file1.replace('_2.fastq.gz', '')) {
                sampleFiles.push(file2.replace('_1.fastq.gz', ''));
            } else if (file1.replace('_R1.fastq.gz', '') === file2.replace('_R2.fastq.gz', '')) {
                const sampleName = file1.replace('_R1.fastq.gz', '');
                sampleFiles.push(sampleName);
                const oldPath1 = path.join(`uploads/${req.uploadFolder}/files`, file1);
                const oldPath2 = path.join(`uploads/${req.uploadFolder}/files`, file2);
                const newPath1 = path.join(`uploads/${req.uploadFolder}/files`, `${sampleName}_1.fastq.gz`);
                const newPath2 = path.join(`uploads/${req.uploadFolder}/files`, `${sampleName}_2.fastq.gz`);

                fs.renameSync(oldPath1, newPath1);
                fs.renameSync(oldPath2, newPath2);
            } else if (file2.replace('_R1.fastq.gz', '') === file1.replace('_R2.fastq.gz', '')) {
                const sampleName = file2.replace('_R1.fastq.gz', '');
                sampleFiles.push(sampleName);
                const oldPath1 = path.join(`uploads/${req.uploadFolder}/files`, file1);
                const oldPath2 = path.join(`uploads/${req.uploadFolder}/files`, file2);
                const newPath1 = path.join(`uploads/${req.uploadFolder}/files`, `${sampleName}_1.fastq.gz`);
                const newPath2 = path.join(`uploads/${req.uploadFolder}/files`, `${sampleName}_2.fastq.gz`);

                fs.renameSync(oldPath1, newPath1);
                fs.renameSync(oldPath2, newPath2);
            } else {
                return res.status(400).send('File pairs do not match.');
            }
        }

        const samplesList = sampleFiles.map(sample => `    "${sample}"`).join(' \\\n');
        const GENOMEIDX1 = indexBasePath;
        const GENOMEIDX = genomeFilePath;

        const basedir = path.resolve(__dirname, `uploads/${req.uploadFolder}`);
        const progressFile = path.join(basedir, 'progress.txt');


        const pipelineScriptContent = `#!/bin/bash

        set -e
        basedir="${basedir}"
        samples=(
        ${samplesList}
        )
        GENOMEIDX1=${GENOMEIDX1}
        GENOMEIDX=${GENOMEIDX}
        progress_file="${basedir}/progress.txt"
        touch "\${progress_file}" && chmod 666 "\${progress_file}"
        echo "Pipeline started" > \${progress_file}

        for sample_name in "\${samples[@]}"; do
            sleep 2
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step I: Quality Control and Preprocessing" >> \${progress_file}
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step 1.1: FastQC Quality Control Report for \${sample_name}" >> \${progress_file}
            fastqc -o "\${basedir}/fastqc_output/" "\${basedir}/files/\${sample_name}_1.fastq.gz" "\${basedir}/files/\${sample_name}_2.fastq.gz" 
            sleep 2
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-1.2: Trimmed Quality Control for \${sample_name}" >> \${progress_file}
            fastp -i "\${basedir}/files/\${sample_name}_1.fastq.gz" -o "\${basedir}/\${sample_name}_P1.fastq" \
            -I "\${basedir}/files/\${sample_name}_2.fastq.gz" -O "\${basedir}/\${sample_name}_P2.fastq" \
            --thread 4 -h "\${basedir}/fastp-\${sample_name}.html" 
            sleep 2
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step II: Read Alignment for \${sample_name}" >> \${progress_file}
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-2.1: Read Alignment for \${sample_name}" >> \${progress_file}
            bowtie2 -p 64 -x "\${GENOMEIDX1}" -1 "\${basedir}/\${sample_name}_P1.fastq" -2 "\${basedir}/\${sample_name}_P2.fastq" -S "\${basedir}/\${sample_name}.sam" 
            sleep 2
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-III: Coverage Analysis for \${sample_name}" >> \${progress_file}
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step 3.1 Conversion Of Sam To BAM File for \${sample_name}" >> \${progress_file}
            samtools view -b "\${basedir}/\${sample_name}.sam" -o "\${basedir}/\${sample_name}.bam" 
            sleep 2
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-3.2: Alignment Metrics for \${sample_name}" >> \${progress_file}
            samtools flagstat "\${basedir}/\${sample_name}.bam" > "\${basedir}/\${sample_name}.flagstat.txt" 
            sleep 2
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-3.3: Conversion of BAM To Sorted BAM for \${sample_name}" >> \${progress_file}
            samtools sort  "\${basedir}/\${sample_name}.bam" -o "\${basedir}/\${sample_name}.sorted.bam" 
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-3.4: Removing duplicate reads from Sorted Bam Files for \${sample_name}" >> \${progress_file} 
            samtools rmdup -S "\${basedir}/\${sample_name}.sorted.bam" "\${basedir}/\${sample_name}.duprem.bam" 
            sleep 2
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-3.5: Deriving Low Coverage Bed File for \${sample_name}" >> \${progress_file} 
            samtools depth "\${basedir}/\${sample_name}.duprem.bam" -aa | awk '$3 < 5 {print $1"\t"$2"\t"$3}' > "\${basedir}/coverage_\${sample_name}.txt" 
            sleep 2
            input_bam="\${basedir}/coverage_\${sample_name}.txt" 
            output_bed="\${basedir}/\${sample_name}.bed" 
            sleep 2
            
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-3.6: Extracting start end coordinates of missing read segments for \${sample_name}" >> \${progress_file}
            python3 "pipelines/convert_txt_to_bed.py" "\${sample_name}" "\${input_bam}" "\${output_bed}" 
            echo "Created \${output_bed} from \${input_txt}" >> \${progress_file}
            sleep 2
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-3.7: Performing N-masking for \${sample_name}" >> \${progress_file}
            bedtools maskfasta -fi "\${GENOMEIDX}" -bed "\${output_bed}" -mc N -fo "\${basedir}/\${sample_name}_masked.fasta"
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step IV: Generation of VCF, VCF Index and Viral Genome for \${sample_name}" >> \${progress_file} 
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-4: Generation of VCF for \${sample_name}" >> \${progress_file}
            bcftools mpileup -f "\${basedir}/\${sample_name}_masked.fasta" "\${basedir}/\${sample_name}.duprem.bam" | bcftools call -cv --ploidy 1 -Oz -o "\${basedir}/\${sample_name}.vcf.gz" 
            sleep 2
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-4.1: Generation of VCF Index for \${sample_name}" >> \${progress_file}
            bcftools index "\${basedir}/\${sample_name}.vcf.gz" 2>&1 >> \${progress_file}
            sleep 2
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Step-4.2: Generation of viral genome fasta for \${sample_name}" >> \${progress_file}
            cat "\${basedir}/\${sample_name}_masked.fasta" | bcftools consensus "\${basedir}/\${sample_name}.vcf.gz" > "\${basedir}/\${sample_name}_genome.fa" 
            sleep 2
            echo "$(date '+%Y-%m-%d %H:%M:%S') - complete!" >> \${progress_file} 
        done
        echo "\$(date '+%Y-%m-%d %H:%M:%S') - Step V: MultiQC Quality Control" >> \${progress_file}
        multiqc "\${basedir}/fastqc_output/" "\${basedir}/" -o "\${basedir}/multiqc_output/" 
                `;

        const scriptPath = path.join(basedir, 'pipeline.sh');
        fs.writeFileSync(scriptPath, pipelineScriptContent);
        fs.chmodSync(scriptPath, '755');

        console.log(`Pipeline script written to ${scriptPath}`);
        sampleFiles.forEach(sample => {
            console.log(`Expecting sample files: ${path.join('uploads', `${sample}_1.fastq.gz`)}, ${path.join('uploads', `${sample}_2.fastq.gz`)}`);
        });

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

        res.json({ message: 'Pipeline execution started.', uploadFolder: req.uploadFolder });
    });
});
// Route to stream progress updates
app.get('/progress/:folder', (req, res) => {
    res.setHeader('Content-Type', 'text/event-stream');
    res.setHeader('Cache-Control', 'no-cache');
    res.setHeader('Connection', 'keep-alive');
    res.flushHeaders();

    const basedir = path.resolve(__dirname, 'uploads/', req.params.folder);
    console.log(basedir)
    const progressFile = path.join(basedir, 'progress.txt');
    let fileOffset = 0;

    const sendProgress = () => {
        if (fs.existsSync(progressFile)) {
            fs.readFile(progressFile, 'utf-8', (err, data) => {
                if (err) {
                    console.error(`Error reading progress file: ${err}`);
                    res.write(`data: Error reading progress file\n\n`);
                    return;
                }

                if (data.length > fileOffset) {
                    const newContent = data.slice(fileOffset);
                    fileOffset = data.length;
                    const newLines = newContent.split('\n').filter(line => line.trim() !== '');
                    newLines.forEach(line => {
                        res.write(`data: ${line}\n\n`);
                    });
                }
            });
        } else {
            res.write('data: Waiting for progress file to be created...\n\n');
        }
    };

    const intervalId = setInterval(sendProgress, 1000);

    req.on('close', () => {
        clearInterval(intervalId);
        res.end();
    });
});


// Route to download the progress log
app.get('/download-progress', (req, res) => {
    const basedir = path.resolve(__dirname, `uploads/${req.params.folder}`);
    const progressFile = path.join(basedir, 'progress.txt');
    if (fs.existsSync(progressFile)) {
        res.download(progressFile, 'progress_log.txt', (err) => {
            if (err) {
                console.error(`Error downloading the file: ${err}`);
                res.status(500).send('Error downloading the file.');
            }
        });
    } else {
        res.status(404).send('Progress log not found.');
    }
});


const UPLOAD_FOLDER = path.join(__dirname, 'uploads');



app.get('/files', (req, res) => {
    fs.readdir(UPLOAD_FOLDER, (err, files) => {
        if (err) {
            return res.status(500).send('Unable to scan files');
        }
        res.json(files);
    });
});

app.get('/uploads/:filename', (req, res) => {
    const file = path.join(UPLOAD_FOLDER, req.params.filename);
    fs.readdir(file, (err, files) => {
        if (err) {
            return res.status(500).send('Unable to scan files');
        }
        res.json(files);
    });
});




app.get('/uploads/:filename/:subfolder', (req, res) => {
    const file = path.join(UPLOAD_FOLDER, req.params.filename, req.params.subfolder);
    fs.readdir(file, (err, files) => {
        if (err) {
            return res.status(500).send('Unable to scan files');
        }
        res.json(files);
    });
});

app.get('/list-directories', (req, res) => {
    const uploadsDir = path.join(__dirname, 'uploads');
    fs.readdir(uploadsDir, (err, files) => {
        if (err) {
            return res.status(500).json({ error: 'Error reading uploads directory' });
        }
        // Filter out non-directory entries and sort by name (which is the timestamp)
        const directories = files.filter(file =>
            fs.statSync(path.join(uploadsDir, file)).isDirectory()
        ).sort((a, b) => b.localeCompare(a)); // Sort in descending order

        res.json(directories);
    });
});

// Update the /retrieve route to serve the HTML file
app.get('/retrieve', (req, res) => {
    res.sendFile(path.join(__dirname, 'views', 'retrieve.html'));
});

// Update the route to serve files
app.get('/download/:folder/:file', (req, res) => {
    const filePath = path.join(__dirname, 'uploads', req.params.folder, req.params.file);
    res.download(filePath, (err) => {
        if (err) {
            res.status(404).send('File not found');
        }
    });
});
// Start the server
app.listen(PORT, () => {
    console.log(`Server is running on http://localhost:${PORT}`);
});
