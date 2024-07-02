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
    return `#!/bin/bash
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
        echo "Step I: Quality Control and Preprocessing" >> \${progress_file}
        echo "Step 1.1: FastQC Quality Control Report for \${sample_name}" >> \${progress_file}
        fastqc -o "\${basedir}/fastqc_output/" "\${basedir}/\${sample_name}_1.fastq.gz" "\${basedir}/\${sample_name}_2.fastq.gz"
        sleep 2
        echo "Step-1.2: Trimmed Quality Control for \${sample_name}" >> \${progress_file}
        fastp -i "\${basedir}/\${sample_name}_1.fastq.gz" -o "\${basedir}/\${sample_name}_P1.fastq" \
              -I "\${basedir}/\${sample_name}_2.fastq.gz" -O "\${basedir}/\${sample_name}_P2.fastq" \
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
    echo "Step 1.3: MultiQC Quality Control" >> \${progress_file}
        multiqc "\${basedir}/fastqc_output/" "\${basedir}/" -o "\${basedir}/multiqc_output/"
    `;

}
// Route to handle file uploads and trigger the pipeline
app.post('/upload', upload.fields([
    { name: 'files', maxCount: 100 },
    { name: 'refgenome', maxCount: 1 }
]), (req, res) => {
    try {
        const files = req.files['files'].sort();
        const refGenome = req.files['refgenome'][0];
        validateFiles(files, refGenome  );

        const date = new Date().toISOString().replace(/:/g, '-');
        const basedir = path.join(UPLOAD_FOLDER, date);
        
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

            res.send('Pipeline execution started.');
        });
    } catch (error) {
        res.status(400).send(error.message);
    }
});

// ... (keep the existing routes for progress, download-progress, and file listing)

// Start the server
app.listen(PORT, () => {
    console.log(`Server is running on http://localhost:${PORT}`);
});