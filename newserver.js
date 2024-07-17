const middleware = require('./utils/middleware');
const bodyParser = require('body-parser');
const multer = require('multer');
const express = require('express');
const bcrypt = require('bcryptjs')
const { exec } = require('child_process');
const fs = require('fs');
const app = express();
const config = require('./utils/config');
const mongoose = require('mongoose');
const User = require('./models/user.model');
const cors = require('cors');
const path = require('path');
// const RedisStore = require('connect-redis')(session)
// const passport = require('passport-local').Strategy
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
app.use(middleware.requestLogger)
app.use(express.static('public'));

app.get('/', (req, res) => {
    res.sendFile(__dirname + '/public/landing.html');
});

// Connect to MongoDB
mongoose.connect('mongodb+srv://shrinidhivasant:shri123@cluster0.qpjxkfz.mongodb.net/?retryWrites=true&w=majority&appName=Cluster0', {
})
.then(() => console.log('Connected to MongoDB'))
.catch(err => console.error('Error connecting to MongoDB:', err));




app.use(express.urlencoded({ extended: true }));

// login and register

// app.post('/register', asyncHandler(async (req, res) => {
//     const { email, name, phone, password, trash } = req.body;

//     if (!name || !email || !password) {
//         res.status(400);
//         throw new Error('Please add all fields');
//     }

//     const userExists = await User.findOne({ email });
//     if (userExists) {
//         res.status(400);
//         throw new Error('User already exists');
//     }

//     const salt = await bcrypt.genSalt(10);
//     const hashedPassword = await bcrypt.hash(password, salt);

//     const user = await User.create({
//         name,
//         email,
//         password: hashedPassword,
//     });

//     if (user) {
//         res.status(201).json({
//             _id: user.id,
//             name: user.name,
//             email: user.email,
//             token: generateToken(user._id),
//         });
//     } else {
//         res.status(400);
//         throw new Error('Invalid user data');
//     }
// }));

// // Login user
// app.post('/login', asyncHandler(async (req, res) => {
//     const { email, password } = req.body;

//     const user = await User.findOne({ email });

//     if (user && (await bcrypt.compare(password, user.password))) {
//         res.json({
//             _id: user.id,
//             name: user.name,
//             email: user.email,
//             token: generateToken(user._id),
//         });
//     } else {
//         res.status(400);
//         throw new Error('Invalid credentials');
//     }
// }));



// // Route to handle form submission
// app.post('/submit-form', async (req, res) => {
//     try {
//         // Create a new document with the form data
        
//         const formData = new FormData({
//             id: (req.body.name.slice(0,4) + req.body.email.slice(0,4) + req.body.mobileno.slice(0,4) + req.body.pwd.slice(0,4)),
//             email: req.body.email,
//             name: req.body.name,
//             mobileNo: req.body.mobileno,
//             password: req.body.pwd
//         });
//         console.log(JSON.stringify(formData));
//         // Save the document to the database
//         await formData.save();

//         res.send('Form data saved to MongoDB!');
//     } catch (err) {
//         console.error('Error saving form data to MongoDB:', err);
//         res.status(500).send('Internal server error');
//     }
// });



const date = new Date().toISOString().replace(/:/g, '-');


const uploadPath = path.join(__dirname, `uploads/${date}/files`);
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


app.post('/upload', upload.array('genomeFiles', 20), (req, res) => {
    const filePaths = req.files.map(file => file.path);
    exec(`bash pipelines/fastqc.sh ${filePaths.join(' ')} ${date}`, (error, stdout, stderr) => {
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

function fetchReports() {
    var xhr = new XMLHttpRequest();
    xhr.open('GET', '/multiqc_reports', true);
    xhr.onload = function () {
        if (xhr.status == 200) {
            var directories = JSON.parse(xhr.responseText);

            // Sort directories by date
            directories.sort(function (a, b) {
                return new Date(b.replace(/-/g, ':')) - new Date(a.replace(/-/g, ':'));
            });

            // Generate the list of links
            var reportListHTML = '<h2>Previous Reports</h2><ul>';
            directories.forEach(function (dir) {
                reportListHTML += `<li><a href="/multiqc_reports/${dir}/multiqc_report.html" target="_blank">${dir}</a></li>`;
            });
            reportListHTML += '</ul>';

            document.getElementById('reportList').innerHTML = reportListHTML;
        } else {
            document.getElementById('reportList').innerHTML = 'Error occurred while fetching directory list';
        }
    };
    xhr.send();
}


app.get('/multiqc_reports', (req, res) => {
    const reportsDir = path.join(__dirname, 'multiqc_reports');
    fs.readdir(reportsDir, (err, files) => {
        if (err) {
            res.status(500).send('Error occurred while fetching reports');
            return;
        }
        res.json(files.filter(file => fs.statSync(path.join(reportsDir, file)).isDirectory()));
    });
});

app.get('/multiqc_report', (req, res) => {
    // Send the MultiQC report as response
    const reportPath = path.join(__dirname, `multiqc_reports/${date}`, 'multiqc_report.html');
    res.sendFile(reportPath);
});

app.listen(config.port, () => {
    console.log(`Server is running on port ${config.port}`);
});