const express = require('express');
const bodyParser = require('body-parser');
const nodemailer = require('nodemailer');

const app = express();
const port = 3000;

app.use(bodyParser.json());
app.use(bodyParser.urlencoded({ extended: true }));

async function sendEmail(to, subject, text) {
    //Create reusable transporter object using the default SMTP transport
    let transporter = nodemailer.createTransport({
        host: 'smtp.gmail.com',
        port: 465,
        secure: true,
        auth: {
            type: 'OAuth2',
            user: process.env.EMAIL,
            clientpassword: process.env.Password
        }
    });
    
    const user_email = 'user-email@gmail.com';
    // Setup email data with unicode symbols
    let mailOptions = {
        from: 'viral-genome-assembly_AIIMSOMICSTOOLBOX@gmail.com', //sender's address (server side)
        to: `${user_email}`, //user's email address (derived from database)    
        subject: 'Analysis Complete',
        text: 'For the following samples submitted, viral genomes are assembled. Click on the link below to download the sample genomes',
        html: `<a href="${link}">Download</a>` //Add href link to download webpage
    }

    //Send email
    transporter.sendEmail(mailOptions, (error, info) => {
       if(error) {
        return console.log(error);
       }
       console.log('Email sent: ', info.response);
    });
}

//Root route
app.get('/', (req, res) => {
  res.send('Automatic Email Response System');
});

// Endpoint to handle email requests
app.post('/', (req, res) => {
    const{to, subject, text} = req.body;

  //Validate request body 
  if (!to ||!subject ||!text) {
    return res.status(400).json({ error: 'Missing required fields' });
    }

   try {
    //Send email
    sendEmail(to, subject, text);
    res.status(200).send('Email successfully sent');
   } catch(error) {
    res.status(500).send('Error sending email: ' + error.message);
   }
});

app.listen(port, () => {
console.log(`Server is running on port ${port}`);
});

//Function to send users email when analysis done
