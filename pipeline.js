const express = require('express');
const multer = require('multer');
const { exec } = require('child_process');
const path = require('path');
const fs = require('fs');

const app = express();
const PORT = process.env.PORT || 3000;

// Set up multer for file uploads
const storage = multer.diskStorage({
    destination: (req, file, cb) => {
        const uploadPath = 'uploads/';
        if (!fs.existsSync(uploadPath)) {
            fs.mkdirSync(uploadPath, { recursive: true });
        }
        cb(null, uploadPath);
    },
    filename: (req, file, cb) => {
        cb(null, file.originalname);
    }
});

const upload = multer({ storage });

// Middleware to serve static files
app.use(express.static(path.join(__dirname, 'public')));

// Route to handle file uploads and trigger the pipeline
app.post('/upload', upload.array('files'), (req, res) => {
    const files = req.files;
    if (files.length % 2 !== 0) {
        return res.status(400).send('Please upload files in pairs.');
    }

    const basedir = path.resolve(__dirname, 'uploads');
    const sampleFiles = [];

    for (let i = 0; i < files.length; i += 2) {
        const file1 = files[i].filename;
        const file2 = files[i + 1].filename;
        if (file1.replace('_1.fastq.gz', '') === file2.replace('_2.fastq.gz', '')) {
            sampleFiles.push(file1.replace('_1.fastq.gz', ''));
        } else {
            return res.status(400).send('File pairs do not match.');
        }
    }

    const samplesList = sampleFiles.map(sample => `    "${sample}"`).join(' \\\n');
    let GENOMEIDX1, GENOMEIDX;

    // Set the variables based on the genome type
    if (req.body.genome === "hev") {
        GENOMEIDX1 = "/home/aryan/omics_website/utils/hev/hev_genome";
        GENOMEIDX = "/home/aryan/omics_website/utils/hev/NC_001434-HEV.fa";
    } else if (req.body.genome === "covid") {
        GENOMEIDX1 = "/home/aryan/omics_website/utils/covid/sars_cov_2";
        GENOMEIDX = "/home/aryan/omics_website/utils/covid/NC_045512.2.fasta";
    } else {
        return res.status(400).send('Invalid genome type specified.');
    }

    const pipelineScriptContent = `#!/bin/bash
    
    set -e
    basedir="${basedir}"
    samples=(
    ${samplesList}
    )
    GENOMEIDX1=${GENOMEIDX1}
    GENOMEIDX=${GENOMEIDX}
    progress_file="${basedir}/progress.txt"

    echo "Pipeline started" > \${progress_file}

    for sample_name in "\${samples[@]}"; do
        sleep 2
        echo "Step-1.0: FastQC Quality Control Report for \${sample_name}" >> \${progress_file}
        fastqc -o "\${basedir}/fastqc_output/" "\${basedir}/\${sample_name}_1.fastq.gz" "\${basedir}/\${sample_name}_2.fastq.gz"
        sleep 2
        echo "Step-1.1: Fastp Quality Control for \${sample_name}" >> \${progress_file}
        fastp -i "\${basedir}/\${sample_name}_1.fastq.gz" -o "\${basedir}/\${sample_name}_P1.fastq" \
              -I "\${basedir}/\${sample_name}_2.fastq.gz" -O "\${basedir}/\${sample_name}_P2.fastq" \
              --thread 4 -h "\${basedir}/fastp-\${sample_name}.html" 2> "\${basedir}/fastp-\${sample_name}.log"
        sleep 2
        echo "Step-2: Read Alignment for \${sample_name}" >> \${progress_file}
        bowtie2 -p 64 -x "\${GENOMEIDX1}" -1 "\${basedir}/\${sample_name}_P1.fastq" -2 "\${basedir}/\${sample_name}_P2.fastq" -S "\${basedir}/\${sample_name}.sam"
        sleep 2
        echo "Step-3: Conversion Of Sam To BAM File for \${sample_name}" >> \${progress_file}
        samtools view -b "\${basedir}/\${sample_name}.sam" -o "\${basedir}/\${sample_name}.bam"
        sleep 2
        echo "Step-4: Alignment Metrics for \${sample_name}" >> \${progress_file}
        samtools flagstat "\${basedir}/\${sample_name}.bam" > "\${basedir}/\${sample_name}.flagstat.txt"
        sleep 2
        echo "Step-5: Conversion of BAM To Sorted BAM for \${sample_name}" >> \${progress_file}
        samtools sort "\${basedir}/\${sample_name}.bam" -o "\${basedir}/\${sample_name}.sorted.bam"
        sleep 2
        echo "Step-6: Deriving Low Coverage Bed File for \${sample_name}" >> \${progress_file}
        samtools depth "\${basedir}/\${sample_name}.sorted.bam" | awk '$3 < 5 {print $1"\t"$2"\t"$3}' > "\${basedir}/coverage_\${sample_name}.txt"
        sleep 2
        input_bam="\${basedir}/coverage_\${sample_name}.txt"
        output_bed="\${basedir}/\${sample_name}.bed"
        sleep 2
        echo "Step-7: Extracting start end coordinates of missing read segments for \${sample_name}" >> \${progress_file}
        python3 "pipelines/convert_bam_to_bed.py" "\${sample_name}" "\${input_bam}" "\${output_bed}"
        sleep 2
        echo "Step-8: Performing N-masking for \${sample_name}" >> \${progress_file}
        bedtools maskfasta -fi "\${GENOMEIDX}" -bed "\${output_bed}" -mc N -fo "\${basedir}/\${sample_name}_masked.fasta"
        sleep 2
        echo "Step-9: Removing duplicate reads from Sorted Bam Files for \${sample_name}" >> \${progress_file}
        samtools rmdup -S "\${basedir}/\${sample_name}.sorted.bam" "\${basedir}/\${sample_name}.duprem.bam"
        sleep 2
        echo "Step-10: Generation of VCF for \${sample_name}" >> \${progress_file}
        bcftools mpileup -f "\${GENOMEIDX}" "\${basedir}/\${sample_name}.duprem.bam" | bcftools call -cv --ploidy 1 -Oz -o "\${basedir}/\${sample_name}.vcf.gz"
        sleep 2
        echo "Step-11: Generation of VCF Index for \${sample_name}" >> \${progress_file}
        bcftools index "\${basedir}/\${sample_name}.vcf.gz"
        sleep 2
        echo "Step-12: Generation of viral genome fasta for \${sample_name}" >> \${progress_file}
        cat "\${basedir}/\${sample_name}_masked.fasta" | bcftools consensus "\${basedir}/\${sample_name}.vcf.gz" > "\${basedir}/\${sample_name}_genome.fa"
        sleep 2
        echo "complete!" >> \${progress_file}
    done
    `;
    
    const scriptPath = path.join(basedir, 'pipeline.sh');
    fs.writeFileSync(scriptPath, pipelineScriptContent);
    fs.chmodSync(scriptPath, '755');
    
    // Debug logs
    console.log(`Pipeline script written to ${scriptPath}`);
    sampleFiles.forEach(sample => {
        console.log(`Expecting sample files: ${path.join(basedir, `${sample}_1.fastq.gz`)}, ${path.join(basedir, `${sample}_2.fastq.gz`)}`);
    });
    
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
    
    res.send('Pipeline execution started.');
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

// Start the server
app.listen(PORT, () => {
    console.log(`Server is running on http://localhost:${PORT}`);
});
