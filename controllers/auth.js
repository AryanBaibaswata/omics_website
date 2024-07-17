const express = require('express');
const jwt = require('jsonwebtoken');
const { User } = require('../models/models');
const authMiddleware = require('../middleware/auth');

const router = express.Router();

// User registration
router.post('/register', async (req, res) => {
    try {
        const user = new User(req.body);
        await user.save();
        const token = jwt.sign({ _id: user._id.toString() }, process.env.JWT_SECRET);
        res.status(201).send({ user, token });
    } catch (error) {
        res.status(400).send(error);
    }
});

// User login
router.post('/login', async (req, res) => {
    try {
        const user = await User.findOne({ email: req.body.email });
        if (!user || !(await user.comparePassword(req.body.password))) {
            return res.status(401).send({ error: 'Login failed!' });
        }
        const token = jwt.sign({ _id: user._id.toString() }, process.env.JWT_SECRET);
        res.send({ user, token });
    } catch (error) {
        res.status(400).send(error);
    }
});

// User logout
router.post('/logout', authMiddleware, async (req, res) => {
    try {
        req.user.tokens = req.user.tokens.filter((token) => token.token !== req.token);
        await req.user.save();
        res.send();
    } catch (error) {
        res.status(500).send();
    }
});

module.exports = router;