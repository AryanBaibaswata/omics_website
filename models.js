const mongoose = require('mongoose');
const { default: reportSchema } = require('./models/reports');
require('dotenv').config();
const url = process.env.MONGODB_URI;
mongoose.set('strictQuery', false);

mongoose.connect(url);

const userSchema = new mongoose.Schema(
    {
        username: {
            type: String,
            required: true,
        },
        password: {
            type: String,
            required: true
        },
        phone: {
            type:String,
        },
        reports,
        pipes,
        issues,
    }
)