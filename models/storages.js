const multer = require('multer');

const storage = (path) => multer.diskStorage({
    destination: (req, file, cb) => {
        const uploadPath = path.join(__dirname, path);
        if (!fs.existsSync(uploadPath)) {
            fs.mkdirSync(uploadPath, { recursive: true });
        }
        cb(null, uploadPath);
    },
    filename: (req, file, cb) => {
        cb(null, file.originalname);
    }
});


module.exports = storage;
