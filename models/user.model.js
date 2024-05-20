const mongoose = require('mongoose')

const userSchema = new mongoose.Schema({
    id: {
        type: String,
        required: true,
    },
    email: {
        type: String,
        required: true
    },
    name: {
        type: String,
        required: true
    },
    mobileNo: {
        type: String,
        required: true
    },
    password: {
        type: String,
        required: true
    },
   });


// userSchema.set('toJSON', {
//     transform: (document, returnedObject) => {

//     }

// })

// Create a model based on the schema
const userModel = mongoose.model('UserData', userSchema);

module.exports = userModel;