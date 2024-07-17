const mongoose = require('mongoose');
const bcrypt = require('bcryptjs');

const userSchema = new mongoose.Schema({
    username: {
        type: String,
        required: true,
        unique: true
    },
    email: {
        type: String,
        required: true,
        unique: true
    },
    password: {
        type: String,
        required: true
    },
    phone: {
        type: String
    },
    reports: [{ type: mongoose.Schema.Types.ObjectId, ref: 'Report' }],
    pipes: [{ type: mongoose.Schema.Types.ObjectId, ref: 'Pipe' }],
    issues: [{ type: mongoose.Schema.Types.ObjectId, ref: 'Issue' }]
});

userSchema.pre('save', async function(next) {
    if (this.isModified('password')) {
        this.password = await bcrypt.hash(this.password, 10);
    }
    next();
});

userSchema.methods.comparePassword = function(candidatePassword) {
    return bcrypt.compare(candidatePassword, this.password);
};

const reportSchema = new mongoose.Schema({
    sample_name: {
        type: [String],
        required: true
    },
    date: {
        type: Date,
        default: Date.now
    },
    flocation: {
        type: String,
        required: true
    },
    sample_type: {
        type: String,
        required: true
    },
    user: {
        type: mongoose.Schema.Types.ObjectId,
        ref: 'User',
        required: true
    },
    status: {
        type: String,
        enum: ['uploaded', 'processing', 'completed', 'failed'],
        default: 'uploaded'
    },
    file_paths: [{
        type: String,
        required: true
    }]
});

const pipeSchema = new mongoose.Schema({
    pipe_params: {
        type: [String],
        required: true
    },
    sample: {
        type: mongoose.Schema.Types.ObjectId,
        ref: 'Report',
        required: true
    },
    user: {
        type: mongoose.Schema.Types.ObjectId,
        ref: 'User',
        required: true
    },
    date: {
        type: Date,
        default: Date.now
    },
    flocation: {
        type: String,
        required: true
    },
    status: {
        type: String,
        enum: ['queued', 'running', 'completed', 'failed'],
        default: 'queued'
    },
    progress: {
        type: Number,
        default: 0,
        min: 0,
        max: 100
    },
    results: {
        type: mongoose.Schema.Types.Mixed
    }
});

const uploadSchema = new mongoose.Schema({
    user: {
        type: mongoose.Schema.Types.ObjectId,
        ref: 'User',
        required: true
    },
    report: {
        type: mongoose.Schema.Types.ObjectId,
        ref: 'Report',
        required: true
    },
    file_name: {
        type: String,
        required: true
    },
    file_size: {
        type: Number,
        required: true
    },
    upload_status: {
        type: String,
        enum: ['pending', 'in_progress', 'completed', 'failed'],
        default: 'pending'
    },
    upload_progress: {
        type: Number,
        default: 0,
        min: 0,
        max: 100
    },
    upload_date: {
        type: Date,
        default: Date.now
    }
});

const User = mongoose.model('User', userSchema);
const Report = mongoose.model('Report', reportSchema);
const Pipe = mongoose.model('Pipe', pipeSchema);
const Upload = mongoose.model('Upload', uploadSchema);

module.exports = { User, Report, Pipe, Upload };