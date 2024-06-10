const mongoose = require('mongoose')

const userSchema = new mongoose.Schema({
    id: {
        type: mongoose.Schema.Types.ObjectId
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
        // required: true
    },
    pwd: {
        type: String,
        required: true
    },
    reports: {
        type: Object
    }, 
    pipes:  {
        type:Object
    }
   });


// userSchema.set('toJSON', {
//     transform: (document, returnedObject) => {

//     }

// })

// Create a model based on the schema
const userModel = mongoose.model('UserData', userSchema);

module.exports = userModel;