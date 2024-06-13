const express = require('express');
const mongoose = require('mongoose');
const app = express();

app.use(express.static('public'));

app.get('/', (req, res) => {
    res.sendFile(__dirname + '/public/login.html');
});

// Connect to MongoDB
mongoose.connect('mongodb+srv://shrinidhivasant:shri123@cluster0.qpjxkfz.mongodb.net/?retryWrites=true&w=majority&appName=Cluster0', {

})
.then(() => console.log('Connected to MongoDB'))
.catch(err => console.error('Error connecting to MongoDB:', err));

// Define a schema for the form data
const formDataSchema = new mongoose.Schema({
    email: String,
    name: String,
    mobileNo: String,
    password: String,
   });

// Create a model based on the schema
const FormData = mongoose.model('FormData', formDataSchema);

app.use(express.urlencoded({ extended: true }));

// Route to handle form submission
app.post('/submit-form', async (req, res) => {
    try {
        // Create a new document with the form data
        const formData = new FormData({
            email: req.body.email,
            name: req.body.name,
            mobileNo: req.body.mobileno,
            password: req.body.pwd
        });

        // Save the document to the database
        await formData.save();

        res.send('Form data saved to MongoDB!');
    } catch (err) {
        console.error('Error saving form data to MongoDB:', err);
        res.status(500).send('Internal server error');
    }
});

const PORT = process.env.PORT || 3000;
app.listen(PORT, () => {
    console.log(`Server is running on port ${PORT}`);
});
