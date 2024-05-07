const express = require('express');
const multer = require('multer');
const mongoose = require('mongoose');
const bodyParser = require('body-parser');

const app = express();
const port = process.env.PORT || 3000; // Set port for deployment

// Configure Multer for file uploads
const upload = multer({ dest: 'uploads/' }); // Change 'uploads/' if needed

// Connect to MongoDB
mongoose.connect('mongodb://your_mongo_connection_string', {
  useNewUrlParser: true,
  useUnifiedTopology: true
})
  .then(() => console.log('Connected to MongoDB'))
  .catch(err => console.error('Error connecting to MongoDB:', err));

// Define Mongoose schema for genome data
const genomeSchema = new mongoose.Schema({
  name: String,
  fastaData: String, // Store the entire FASTA data as a string for simplicity
  fastqcReport: String // Store the FASTQC report as a string
});

const Genome = mongoose.model('Genome', genomeSchema);

// Configure body-parser for form data access
app.use(bodyParser.urlencoded({ extended: false }));

// Route for uploading genome file
app.post('/upload-genome', upload.single('genomeFile'), async (req, res) => {
  try {
    const { filename, originalname } = req.file;
    const fastaData = await readFileSync(`uploads/${filename}`, 'utf-8'); // Read uploaded FASTA data

    const genome = new Genome({
      name: originalname,
      fastaData,
      fastqcReport: '' // Initialize empty report field
    });

    await genome.save();

    // Call the Bash script to run FASTQC (see below)
    const fastqcReport = await execSync(`bash fastqc_pipeline.sh ${filename}`); // Modify command if needed

    // Update the genome data with the report
    await Genome.updateOne({ _id: genome._id }, { fastqcReport: fastqcReport.toString() });

    res.send('Genome uploaded and processed successfully!');
  } catch (err) {
    console.error('Error uploading or processing genome:', err);
    res.status(500).send('Error processing genome. Please try again.');
  }
});

// Route to serve FASTQC report (assuming the report is stored in the database)
app.get('/report/:genomeId', async (req, res) => {
  try {
    const genome = await Genome.findById(req.params.genomeId);
    if (!genome) {
      return res.status(404).send('Genome not found.');
    }

    res.send(genome.fastqcReport); // Send the report data
  } catch (err) {
    console.error('Error fetching FASTQC report:', err);
    res.status(500).send('Error retrieving report. Please try again.');
  }
});

app.listen(port, () => console.log(`Server listening on port ${port}`));
