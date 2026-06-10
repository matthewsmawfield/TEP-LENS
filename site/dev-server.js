#!/usr/bin/env node
const chokidar = require('chokidar');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const net = require('net');
const { buildStaticSite } = require('./build.js');
class DevServer {
    constructor() {
        this.isBuilding = false;
        this.buildQueue = false;
        this.liveServerProcess = null;
        this.watcherReady = false;
        this.port = 51736; // Unique port for TEP-LENS
    }

    async killProcessOnPort(port) {
        return new Promise((resolve) => {
            const cleanup = spawn('sh', ['-c', `lsof -ti :${port} | xargs kill -9 2>/dev/null || true`], { stdio: 'ignore' });
            cleanup.on('close', () => setTimeout(resolve, 500));
        });
    }

    async startLiveServer() {
        console.log('🚀 Starting live server...');
        if (this.liveServerProcess) this.liveServerProcess.kill('SIGKILL');
        await this.killProcessOnPort(this.port);
        this.liveServerProcess = spawn('npx', [
            'live-server', 'dist', `--port=${this.port}`, '--host=localhost', '--no-browser', '--wait=500'
        ], { stdio: 'pipe', cwd: __dirname });
        this.liveServerProcess.stdout.on('data', d => console.log('📡 ' + d.toString().trim()));
    }

    async build() {
        if (this.isBuilding) { this.buildQueue = true; return; }
        this.isBuilding = true;
        console.log('\n🔄 Rebuilding site...');
        try {
            await buildStaticSite();
            console.log('✅ Build complete (includes markdown generation)');
            if (this.buildQueue) {
                this.buildQueue = false;
                setTimeout(() => { this.isBuilding = false; this.build(); }, 100);
                return;
            }
        } catch (e) { console.error('❌ Build failed:', e.message); }
        this.isBuilding = false;
    }

    async start() {
        process.on('SIGINT', () => {
            if (this.liveServerProcess) this.liveServerProcess.kill('SIGKILL');
            process.exit(0);
        });
        console.log('🎯 TEP-LENS Development Server');
        const distDir = path.join(__dirname, 'dist');
        if (!fs.existsSync(distDir)) fs.mkdirSync(distDir, { recursive: true });
        await this.build();
        await this.startLiveServer();
        
        const watcher = chokidar.watch([
            path.join(__dirname, 'components'),
            path.join(__dirname, 'index.html'),
            path.join(__dirname, 'journal.html'),
            path.join(__dirname, 'manifest.json'),
            path.join(__dirname, 'styles'),
            path.join(__dirname, 'figures')
        ], { ignored: ['dist/**'], persistent: true, ignoreInitial: true });

        watcher.on('change', () => this.build());
        watcher.on('add', () => this.build());
        console.log(`\n🌐 Server running at: http://localhost:${this.port}`);
    }
}

if (require.main === module) { const s = new DevServer(); s.start(); }
module.exports = DevServer;
