# Fusang Web Server - Conda Deployment Guide
# ==========================================

## 1. Prerequisites

- **Conda** (Miniconda or Anaconda) installed
  - Download: https://docs.conda.io/en/latest/miniconda.html
- Minimum recommended: 4 GB RAM, 2 CPU cores

## 2. Quick Start (Conda Environment)

```bash
# Step 1: Navigate to the Fusang project directory
cd /path/to/Fusang-main

# Step 2: Create conda environment from environment.yml
conda env create -f environment.yml

# Step 3: Activate the environment
conda activate fusang-web

# Step 4: Start the web server
fusang-web

# Step 5: Open browser
# → http://localhost:5000
```

## 3. Configuration Options

### Custom Host / Port

```bash
# Start on port 8080, local only
fusang-web --host 127.0.0.1 --port 8080

# Start on all interfaces (accessible from LAN)
fusang-web --host 0.0.0.0 --port 5000

# Enable debug mode (development only)
fusang-web --debug
```

### Environment Variables

| Variable            | Default    | Description                  |
|---------------------|------------|------------------------------|
| `FUSANG_HOST`       | `0.0.0.0`  | Bind address                 |
| `FUSANG_PORT`       | `5000`     | Listen port                  |
| `FUSANG_DEBUG`      | `false`    | Enable Flask debug mode      |
| `FUSANG_MAX_UPLOAD` | `500`      | Max upload size in MB        |
| `FUSANG_THREADS`    | `4`        | Default CPU threads for jobs |

## 4. Production Deployment

### Option A: Gunicorn (Linux/macOS, recommended)

```bash
# Install gunicorn (included in environment.yml)
conda activate fusang-web
gunicorn -w 4 -b 0.0.0.0:5000 fusang_webapp:app

# With logging
gunicorn -w 4 -b 0.0.0.0:5000 \
  --access-logfile /var/log/fusang/access.log \
  --error-logfile /var/log/fusang/error.log \
  fusang_webapp:app
```

### Option B: Waitress (Cross-platform, Windows/Linux/macOS)

```bash
# Install waitress (included in environment.yml)
conda activate fusang-web
waitress-serve --host=0.0.0.0 --port=5000 fusang_webapp:app
```

### Option C: Systemd Service (Linux)

Create `/etc/systemd/system/fusang-web.service`:

```ini
[Unit]
Description=Fusang Web Server
After=network.target

[Service]
Type=simple
User=fusang
Group=fusang
WorkingDirectory=/opt/fusang
Environment="PATH=/opt/conda/envs/fusang-web/bin"
ExecStart=/opt/conda/envs/fusang-web/bin/gunicorn -w 4 -b 0.0.0.0:5000 fusang_webapp:app
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

Then:
```bash
sudo systemctl daemon-reload
sudo systemctl enable fusang-web
sudo systemctl start fusang-web
sudo systemctl status fusang-web
```

## 5. Building the Conda Package

```bash
# Install conda-build
conda install conda-build

# Build the package
cd Fusang-main
conda build conda-recipe/

# Install the built package locally
conda install --use-local fusang-web
```

## 6. Updating

```bash
# Update the conda environment
conda activate fusang-web
conda env update -f environment.yml --prune

# Or update individual packages
conda update --all -n fusang-web
```

## 7. File Structure

```
Fusang-main/
├── conda-recipe/           # Conda build recipe
│   ├── meta.yaml           # Package metadata and dependencies
│   ├── build.sh            # Unix build script
│   └── bld.bat             # Windows build script
├── environment.yml         # Conda environment definition
├── setup.py                # pip installer with entry point
├── MANIFEST.in             # Include non-Python files
├── fusang_webapp.py        # Web server (Flask)
├── fusang_v2.py            # Core Fusang engine
├── fusang_mhl/             # MHL package
├── templates/              # HTML templates
│   ├── index.html          # Upload page
│   ├── result.html         # Job status page
│   └── visualize.html      # Tree visualization (D3.js)
├── static/                 # Static assets + results
│   ├── results/            # Generated Newick trees
│   ├── css/
│   └── js/
└── DEPLOYMENT.md           # This file
```

## 8. Troubleshooting

### Port already in use
```bash
# Find and kill the process using the port
lsof -i :5000
kill -9 <PID>
```

### Conda environment creation fails
```bash
# Create empty env first, then install
conda create -n fusang-web python=3.11
conda activate fusang-web
pip install -e /path/to/Fusang-main
```

### "Module not found" errors
```bash
# Ensure the package is installed in development mode
conda activate fusang-web
pip install -e /path/to/Fusang-main
```

### Permission denied (binding to port < 1024)
```bash
# Use a port >= 1024, or use authbind (Linux)
fusang-web --port 8080
```

## 9. Security Notes

- `--debug` mode exposes detailed error messages — never use in production
- The Flask development server (`fusang-web`) is for development/testing only
- For production, use Gunicorn/Waitress behind Nginx reverse proxy
- Configure firewall to restrict access if needed
- Default maximum upload is 500 MB — adjust in `fusang_webapp.py` line 37
