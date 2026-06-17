#!/bin/bash
# ============================================================
# Fusang Web Server Launcher (Linux/macOS)
# ============================================================
# Automatically activates conda env and starts the web server.
# For production, use: gunicorn -w 4 -b 0.0.0.0:5000 fusang_webapp:app
# ============================================================

set -e

ENV_NAME="fusang-web"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Activate conda environment
if ! conda activate "${ENV_NAME}" 2>/dev/null; then
    echo "[ERROR] Cannot activate conda environment '${ENV_NAME}'."
    echo "Run 'bash setup_conda.sh' first to create the environment."
    exit 1
fi

echo ""
echo "============================================================"
echo "  Fusang Web Server"
echo "  Access at: http://localhost:5000"
echo "============================================================"
echo ""

cd "${SCRIPT_DIR}"

if [ "${1}" = "--prod" ] || [ "${1}" = "-p" ]; then
    # Production mode with Gunicorn
    echo "[INFO] Starting in production mode with Gunicorn (4 workers)..."
    gunicorn -w 4 -b 0.0.0.0:5000 fusang_webapp:app
else
    # Development mode (pass through CLI args)
    python fusang_webapp.py "$@"
fi
