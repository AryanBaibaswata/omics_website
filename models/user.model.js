const userModelSchema = new mongoose.Schema({
    id: String,
    email: String,
    name: String,
    mobileNo: String,
    password: String,
   });


// Create a model based on the schema
const userModel = mongoose.model('UserData', userModelSchema);