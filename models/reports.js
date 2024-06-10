const mongoose = require("mongoose");

const reportSchema = new Schema({
    userID: { type: String, required: true },
    date: { type: Date, default: Date.now, required: true },
    fileID: { type: mongoose.Schema.Types.ObjectId, required: true, ref: 'fs.files' },
  });

const Report = mongoose.model('Report', reportSchema);

export default reportSchema;