const mongoose = require('mongoose');
const { default: reportSchema } = require('./models/reports');
const bcrypt = require('bcryptjs');
require('dotenv').config();
const url = process.env.MONGODB_URI;
mongoose.set('strictQuery', false);

mongoose.connect(url);

const userSchema = new mongoose.Schema(
    {
        user_id: {
            type: mongoose.Schema.Types.ObjectId,
            required: true,
        },
        username: {
            type: String,
            required: true,
        },
        password: {
            type: String,
            required: true,
        },
        phone: {
            type:String,
        },
        reports,
        pipes,
        issues,
    }
)

const reportSchema = new mongoose.Schema(
    {
        sample_id: {
            type: mongoose.Schema.Types.ObjectId,
            required: true,
        },
        sample_name:{
            type: Array,
            required: true,
        },
        date: {
            type: Date,
            required: true,
        },
        flocation: {
            type: String,
            required: true,
        },
        sample_type: {
            type: String,
            required: true,
        },
        user_id: {
            type: mongoose.Schema.Types.ObjectId,
            required: true,
        }
    }
)

const pipeSchema = new mongoose.Schema(
    {
        pipe_id: {
            type: mongoose.Schema.Types.ObjectId,
            required: true,
        },
        pipe_params:{
            type: Array,
            required: true,
        },
        sample_id: {
            type: mongoose.Schema.Types.ObjectId,
            required: true,
        },
        user_id: {
            type: mongoose.Schema.Types.ObjectId,
            required: true,
        },
        date: {
            type: Date,
            required: true,
        },
        flocation: {
            type: String,
            required: true,
        },
        
    }
)

userSchema.pre('save', async function(next) {
    if(!this.isModified('password')) {
        return next();
    }
    const salt = await bcrypt.genSalt(10);
    this.password = await bcrypt.hash(this.password, salt);
    next();
});

userSchema.methods.matchPassword = async function(enteredPassword) {
    return await bcrypt.compare(enteredPassword, this.password);
};




module.exports = {
    User: mongoose.model('User', userSchema),
    Report: mongoose.model('Report', reportSchema),
    Pipe: mongoose.model('Pipe', pipeSchema),
}