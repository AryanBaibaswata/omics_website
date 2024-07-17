const app = require('express');

app.use(express.static(path.join(__dirname, 'public')));
