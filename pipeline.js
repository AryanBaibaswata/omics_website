const express = require('express');
const multer = require('multer');
const { exec } = require('child_process');
const path = require('path');
const fs = require('fs');
const json = require('json');
// // const passport = require('passport');
// // const LocalStrategy = require('passport-local').Strategy;
// // const session = require('express-session');
// // const redis = require("connect-redis")(session);
// // const bcrypt = require('bcrypt');
const app = express();
const PORT = process.env.PORT || 3000;
// const handlebars = require('express-handlebars');
// const handlebars = require('express-handlebars');
// Set up multer for file uploads

const date = new Date().toISOString().replace(/:/g, '-');

const storage = multer.diskStorage({
    destination: (req, file, cb) => {
        const uploadPath = `uploads/${date}/files`;
        if (!fs.existsSync(uploadPath)) {
            fs.mkdirSync(uploadPath, { recursive: true });
            fs.mkdirSync(`uploads/${date}/fastqc_output`, { recursive: true });
        }
        cb(null, uploadPath);
    },
    filename: (req, file, cb) => {
        cb(null, file.originalname);
    }
});


// app.use(session({
//     store: new RedisStore({
//       url: config.redisStore.url
//     }),
//     secret: config.redisStore.secret,
//     resave: false,
//     saveUninitialized: false
//   }))
//   app.use(passport.initialize())
//   app.use(passport.session())


// app.engine('.hbs', handlebars.engine({
//     defaultLayout: 'main',
//     extname: '.hbs',
//     layoutsDir: path.join(__dirname, 'views/layouts')
// }));
const upload = multer({ storage }).fields([
    { name: 'files' }, // Adjust maxCount as needed
    { name: 'refgenome', maxCount: 1 }
]);
// console.log("read part till storage")

// Middleware to serve static files
app.use(express.static(path.join(__dirname, 'public')));
// console.log("read till middleware")
// Route to handle file uploads and trigger the pipeline
app.post('/upload', upload, (req, res) => {
    const files = req.files['files'];
    const genomeFiles = req.files['refgenome'];

    if (!genomeFiles || genomeFiles.length === 0) {
        return res.status(400).send('Please upload a .fasta or .fa file.');
    }

    const genomeFile = genomeFiles[0];
    const genomeFilePath = path.join(`uploads/${date}/files`, genomeFile.filename);
    console.log(genomeFilePath);
    if (files.length % 2 !== 0) {
        return res.status(400).send('Please upload files in pairs.');
    }
    const genomeFilename = path.basename(genomeFile.filename, path.extname(genomeFile.filename));
    console.log(genomeFilename)
    fs.mkdirSync(`uploads/${date}/${genomeFilename}`,)
    const indexBasePath = path.join(`uploads/${date}/${genomeFilename}`, path.basename(genomeFile.filename, path.extname(genomeFile.filename)));
    console.log("indexbasepath:", indexBasePath)
    const bowtieBuildCmd = `bowtie2-build ${genomeFilePath} ${indexBasePath}`;

    exec(bowtieBuildCmd, (error, stdout, stderr) => {
        console.log("began execing")
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
                const oldPath1 = path.join(`uploads/${date}/files`, file1);
                const oldPath2 = path.join(`uploads/${date}/files`, file2);
                const newPath1 = path.join(`uploads/${date}/files`, `${sampleName}_1.fastq.gz`);
                const newPath2 = path.join(`uploads/${date}/files`, `${sampleName}_2.fastq.gz`);

                fs.renameSync(oldPath1, newPath1);
                fs.renameSync(oldPath2, newPath2);
            } else if (file2.replace('_R1.fastq.gz', '') === file1.replace('_R2.fastq.gz', '')) {
                const sampleName = file2.replace('_R1.fastq.gz', '');
                sampleFiles.push(sampleName);
                const oldPath1 = path.join(`uploads/${date}/files`, file1);
                const oldPath2 = path.join(`uploads/${date}/files`, file2);
                const newPath1 = path.join(`uploads/${date}/files`, `${sampleName}_1.fastq.gz`);
                const newPath2 = path.join(`uploads/${date}/files`, `${sampleName}_2.fastq.gz`);

                fs.renameSync(oldPath1, newPath1);
                fs.renameSync(oldPath2, newPath2);
            } else {
                return res.status(400).send('File pairs do not match.');
            }
        }

        const samplesList = sampleFiles.map(sample => `    "${sample}"`).join(' \\\n');
        const GENOMEIDX1 = indexBasePath;
        const GENOMEIDX = genomeFilePath;
        const basedir = path.resolve(__dirname, `uploads/${date}`);


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
            echo "Step I: Quality Control and Preprocessing" >> \${progress_file}
            echo "Step 1.1: FastQC Quality Control Report for \${sample_name}" >> \${progress_file}
            fastqc -o "\${basedir}/fastqc_output/" "\${basedir}/files/\${sample_name}_1.fastq.gz" "\${basedir}/files/\${sample_name}_2.fastq.gz"
            sleep 2
            echo "Step-1.2: Trimmed Quality Control for \${sample_name}" >> \${progress_file}
            fastp -i "\${basedir}/files/\${sample_name}_1.fastq.gz" -o "\${basedir}/\${sample_name}_P1.fastq" \
                  -I "\${basedir}/files/\${sample_name}_2.fastq.gz" -O "\${basedir}/\${sample_name}_P2.fastq" \
                  --thread 4 -h "\${basedir}/fastp-\${sample_name}.html" 2> "\${basedir}/fastp-\${sample_name}.log"
            sleep 2
            echo "Step II: Read Alignment for \${sample_name}" >> \${progress_file}
            echo "Step-2.1: Read Alignment for \${sample_name}" >> \${progress_file}
            bowtie2 -p 64 -x "\${GENOMEIDX1}" -1 "\${basedir}/\${sample_name}_P1.fastq" -2 "\${basedir}/\${sample_name}_P2.fastq" -S "\${basedir}/\${sample_name}.sam"
            sleep 2
            echo "Step-III: Coverage Analysis for \${sample_name}" >> \${progress_file}
            echo "Step 3.1 Conversion Of Sam To BAM File for \${sample_name}" >> \${progress_file}
            samtools view -b "\${basedir}/\${sample_name}.sam" -o "\${basedir}/\${sample_name}.bam"
            sleep 2
            echo "Step-3.2: Alignment Metrics for \${sample_name}" >> \${progress_file}
              samtools flagstat "\${basedir}/\${sample_name}.bam" > "\${basedir}/\${sample_name}.flagstat.txt"
            sleep 2
            echo "Step-3.3: Conversion of BAM To Sorted BAM for \${sample_name}" >> \${progress_file}
            samtools sort "\${basedir}/\${sample_name}.bam" -o "\${basedir}/\${sample_name}.sorted.bam"
            echo "Step-3.4: Removing duplicate reads from Sorted Bam Files for \${sample_name}" >> \${progress_file}
            samtools rmdup -S "\${basedir}/\${sample_name}.sorted.bam" "\${basedir}/\${sample_name}.duprem.bam"
            sleep 2
            echo "Step-3.5: Deriving Low Coverage Bed File for \${sample_name}" >> \${progress_file}
            samtools depth "\${basedir}/\${sample_name}.duprem.bam" | awk '$3 < 5 {print $1"\t"$2"\t"$3}' > "\${basedir}/coverage_\${sample_name}.txt"
            sleep 2
            input_bam="\${basedir}/coverage_\${sample_name}.txt"
            output_bed="\${basedir}/\${sample_name}.bed"
            sleep 2
            echo "Step-3.6: Extracting start end coordinates of missing read segments for \${sample_name}" >> \${progress_file}
            python3 "pipelines/convert_bam_to_bed.py" "\${sample_name}" "\${input_bam}" "\${output_bed}"
            sleep 2
            echo "Step-3.7: Performing N-masking for \${sample_name}" >> \${progress_file}
            bedtools maskfasta -fi "\${GENOMEIDX}" -bed "\${output_bed}" -mc N -fo "\${basedir}/\${sample_name}_masked.fasta"
            echo "Step IV: Generation of VCF, VCF Index and Viral Genome for \${sample_name}" >> \${progress_file}
            echo "Step-4: Generation of VCF for \${sample_name}" >> \${progress_file}
            bcftools mpileup -f "\${GENOMEIDX}" "\${basedir}/\${sample_name}.duprem.bam" | bcftools call -cv --ploidy 1 -Oz -o "\${basedir}/\${sample_name}.vcf.gz"
            sleep 2
            echo "Step-4.1: Generation of VCF Index for \${sample_name}" >> \${progress_file}
            bcftools index "\${basedir}/\${sample_name}.vcf.gz"
            sleep 2
            echo "Step-4.2: Generation of viral genome fasta for \${sample_name}" >> \${progress_file}
            cat "\${basedir}/\${sample_name}_masked.fasta" | bcftools consensus "\${basedir}/\${sample_name}.vcf.gz" > "\${basedir}/\${sample_name}_genome.fa"
            sleep 2
            echo "complete!" >> \${progress_file}
        done
        echo \Step 1.3: MultiQC Quality Control" >> \${progress_file}
            multiqc "\${basedir}/fastqc_output/" "\${basedir}/" -o "\${basedir}/multiqc_output/"
        `;

        const scriptPath = path.join(basedir, 'pipeline.sh');
        fs.writeFileSync(scriptPath, pipelineScriptContent);
        fs.chmodSync(scriptPath, '755');

        // Debug logs
        console.log(`Pipeline script written to ${scriptPath}`);
        sampleFiles.forEach(sample => {
            console.log(`Expecting sample files: ${path.join('uploads', `${sample}_1.fastq.gz`)}, ${path.join('uploads', `${sample}_2.fastq.gz`)}`);
        });

        const progressFile = path.join('uploads', 'progress.txt');
        fs.writeFileSync(progressFile, ''); // Clear previous progress

        const child = exec(`stdbuf -oL -eL bash ${scriptPath}`, { shell: '/bin/bash' });
        child.stdout.pipe(process.stdout);
        child.stderr.pipe(process.stderr);

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


        res.send('Pipeline execution started.');
    });
});
// Route to stream progress updates
app.get('/progress', (req, res) => {
    res.setHeader('Content-Type', 'text/event-stream');
    res.setHeader('Cache-Control', 'no-cache');
    res.setHeader('Connection', 'keep-alive');
    res.flushHeaders(); // flush the headers to establish SSE with client

    const progressFile = path.join(__dirname, 'uploads', 'progress.txt');
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
                    const newLines = newContent.split('\n').filter(line => line.trim() !== '').join('\n');
                    if (newLines) {
                        res.write(`data: ${newLines}\n\n`);
                    }
                }
            });
        } else {
            res.write('data: Waiting for progress file to be created...\n\n');
        }
    };

    const intervalId = setInterval(sendProgress, 500); // Reduced interval to capture updates more frequently

    req.on('close', () => {
        clearInterval(intervalId);
        res.end();
    });
});

// Route to download the progress log
app.get('/download-progress', (req, res) => {
    const progressFile = path.join(__dirname, 'uploads', 'progress.txt');
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
console.log

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
// Start the server
app.listen(PORT, () => {
    console.log(`Server is running on http://localhost:${PORT}`);
});
