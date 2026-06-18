"""
Fusang Web Application
Web interface for Fusang phylogenetic inference tool.

Usage:
    # Direct run
    python fusang_webapp.py [--host 0.0.0.0] [--port 5000] [--debug]

    # After conda/pip install
    fusang-web [--host 0.0.0.0] [--port 5000] [--debug]

    # Production (after conda env activation)
    gunicorn -w 4 -b 0.0.0.0:5000 fusang_webapp:app
"""

import os
import sys
import json
import time
import uuid
import threading
import argparse
import logging
from pathlib import Path

from flask import Flask, render_template, request, jsonify, send_file, url_for, redirect
import subprocess

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('fusang-web')

# ===========================================================
# Configuration
# ===========================================================

BASE_DIR = Path(__file__).parent
UPLOAD_FOLDER = BASE_DIR / "uploads"
RESULT_FOLDER = BASE_DIR / "static" / "results"
TEMPLATE_FOLDER = BASE_DIR / "templates"

# Create directories if not exist
UPLOAD_FOLDER.mkdir(exist_ok=True)
RESULT_FOLDER.mkdir(parents=True, exist_ok=True)

app = Flask(__name__, 
            template_folder=str(TEMPLATE_FOLDER),
            static_folder=str(BASE_DIR / "static"))

app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024  # 500MB max
app.config['UPLOAD_FOLDER'] = str(UPLOAD_FOLDER)
app.config['RESULT_FOLDER'] = str(RESULT_FOLDER)

# Job storage (in-memory, for production use Redis/Database)
jobs = {}

# ===========================================================
# Background Task Processing
# ===========================================================

def run_fusang_job(job_id, input_file, params):
    """
    Run Fusang in background thread.
    Updates job status as it progresses.
    """
    job = jobs[job_id]
    job['status'] = 'running'
    job['start_time'] = time.time()
    
    try:
        # Generate output filename
        output_file = RESULT_FOLDER / f"{job_id}.nwk"
        
        # Build command
        cmd = [
            sys.executable,
            str(BASE_DIR / "fusang_v2.py"),
            "-i", str(input_file),
            "-o", str(output_file),
            "-m", params.get('mode', 'auto'),
            "-d", params.get('distance_method', 'kmer'),
            "-t", params.get('threads', '4'),
            "--tree_method", params.get('tree_method', 'nj'),
            "--max_group", params.get('max_group', '200'),
            "--overlap", params.get('overlap', '0.15'),
        ]
        
        # Add optional parameters
        if params.get('kmer_k'):
            cmd.extend(["--kmer_k", str(params['kmer_k'])])
        if params.get('kmer_gap'):
            cmd.extend(["--kmer_gap", params['kmer_gap']])
        if params.get('simple') == 'true':
            cmd.append("--simple")
            
        job['command'] = ' '.join(cmd)
        job['log'].append(f"[INFO] Starting Fusang...")
        job['log'].append(f"[INFO] Command: {' '.join(cmd)}")
        
        # Run subprocess
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=str(BASE_DIR)
        )
        
        # Capture output
        output_lines = []
        for line in iter(process.stdout.readline, ''):
            if line:
                output_lines.append(line.strip())
                job['log'].append(line.strip())
                # Update progress based on output
                if '[MAIN]' in line or '[kmer]' in line or '[DCM]' in line or '[NNI]' in line:
                    job['progress'] = min(95, job['progress'] + 2)
        
        process.wait()
        
        if process.returncode == 0:
            job['status'] = 'completed'
            job['progress'] = 100
            job['output_file'] = str(output_file)
            job['log'].append(f"[INFO] Job completed successfully!")
            job['log'].append(f"[INFO] Output: {output_file}")
        else:
            job['status'] = 'failed'
            job['error'] = f"Fusang exited with code {process.returncode}"
            job['log'].append(f"[ERROR] Fusang failed with return code {process.returncode}")
            
    except Exception as e:
        job['status'] = 'failed'
        job['error'] = str(e)
        job['log'].append(f"[ERROR] Exception: {str(e)}")
    
    finally:
        job['end_time'] = time.time()
        job['duration'] = round(job['end_time'] - job['start_time'], 2)


# ===========================================================
# Routes
# ===========================================================

@app.route('/')
def index():
    """Home page - upload form."""
    return render_template('index.html')


@app.route('/upload', methods=['POST'])
def upload():
    """Handle file upload and start processing."""
    try:
        # Check file
        if 'fasta_file' not in request.files:
            return jsonify({'error': 'No file uploaded'}), 400
        
        file = request.files['fasta_file']
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        
        # Save uploaded file
        job_id = str(uuid.uuid4())
        filename = f"{job_id}_{file.filename}"
        filepath = UPLOAD_FOLDER / filename
        file.save(str(filepath))
        
        # Get parameters
        params = {
            'mode': request.form.get('mode', 'auto'),
            'distance_method': request.form.get('distance_method', 'kmer'),
            'kmer_k': request.form.get('kmer_k', ''),
            'kmer_gap': request.form.get('kmer_gap', 'none'),
            'tree_method': request.form.get('tree_method', 'nj'),
            'max_group': request.form.get('max_group', '200'),
            'overlap': request.form.get('overlap', '0.15'),
            'threads': request.form.get('threads', '4'),
            'simple': request.form.get('simple', 'false'),
        }
        
        # Clean params
        if params['kmer_k'] == '':
            params['kmer_k'] = None
        if params['kmer_gap'] == 'none':
            params['kmer_gap'] = None
            
        # Create job
        jobs[job_id] = {
            'id': job_id,
            'status': 'pending',
            'progress': 0,
            'params': params,
            'input_file': str(filepath),
            'output_file': None,
            'log': [],
            'error': None,
            'start_time': None,
            'end_time': None,
            'duration': None,
            'filename': file.filename,
        }
        
        # Start background thread
        thread = threading.Thread(
            target=run_fusang_job,
            args=(job_id, filepath, params)
        )
        thread.daemon = True
        thread.start()
        
        return jsonify({
            'job_id': job_id,
            'redirect': url_for('result', job_id=job_id)
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/result/<job_id>')
def result(job_id):
    """Show job result page."""
    job = jobs.get(job_id)
    if not job:
        return "Job not found", 404
    
    return render_template('result.html', job=job)


@app.route('/api/job/<job_id>')
def api_job_status(job_id):
    """API endpoint to get job status (for AJAX polling)."""
    job = jobs.get(job_id)
    if not job:
        return jsonify({'error': 'Job not found'}), 404
    
    return jsonify({
        'id': job['id'],
        'status': job['status'],
        'progress': job['progress'],
        'log': job['log'][-50:],  # Last 50 lines
        'error': job['error'],
        'output_file': job.get('output_file'),
        'duration': job.get('duration'),
    })


@app.route('/download/<job_id>')
def download(job_id):
    """Download result file."""
    job = jobs.get(job_id)
    if not job or not job.get('output_file'):
        return "File not found", 404
    
    output_file = job['output_file']
    if not os.path.exists(output_file):
        return "File not found", 404
    
    return send_file(
        output_file,
        as_attachment=True,
        download_name=f"fusang_tree_{job_id}.nwk",
        mimetype='text/plain'
    )


@app.route('/visualize/<job_id>')
def visualize(job_id):
    """Visualize tree using D3.js."""
    job = jobs.get(job_id)
    if not job or not job.get('output_file'):
        return "Job or result not found", 404
    
    # Read Newick file
    try:
        with open(job['output_file'], 'r') as f:
            newick_str = f.read().strip()
        return render_template('visualize.html', 
                               job_id=job_id, 
                               newick=newick_str)
    except Exception as e:
        return f"Error reading tree file: {str(e)}", 500


@app.route('/jobs')
def list_jobs():
    """List all jobs (for debugging)."""
    return jsonify({
        'jobs': [
            {
                'id': j['id'],
                'status': j['status'],
                'filename': j['filename'],
                'progress': j['progress'],
            }
            for j in jobs.values()
        ]
    })


# ===========================================================
# CLI and Main
# ===========================================================

def parse_args():
    """Parse command-line arguments for the web server."""
    parser = argparse.ArgumentParser(
        description="Fusang Web Server - Scalable Phylogenetic Inference",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  fusang-web                                    # Start on default port 5000
  fusang-web --port 8080                        # Custom port
  fusang-web --host 127.0.0.1 --port 9090       # Local only on port 9090
  fusang-web --debug                            # Enable debug mode
  gunicorn -w 4 -b 0.0.0.0:5000 fusang_webapp:app  # Production (conda env)
        """
    )
    parser.add_argument(
        "--host", type=str, default="0.0.0.0",
        help="Host to bind to (default: 0.0.0.0, all interfaces)"
    )
    parser.add_argument(
        "--port", type=int, default=5000,
        help="Port to listen on (default: 5000)"
    )
    parser.add_argument(
        "--debug", action="store_true", default=False,
        help="Enable Flask debug mode (NOT for production)"
    )
    parser.add_argument(
        "--version", action="version", version="fusang-web 1.0.0"
    )
    return parser.parse_args()


def main():
    """Entry point for fusang-web console script."""
    args = parse_args()

    logger.info("=" * 60)
    logger.info("  Fusang Web Server v1.0.0")
    logger.info("  Host: %s  Port: %d", args.host, args.port)
    logger.info("  Debug: %s", args.debug)
    logger.info("  Access at: http://%s:%d", args.host.replace("0.0.0.0", "localhost"), args.port)
    logger.info("=" * 60)

    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == '__main__':
    main()
