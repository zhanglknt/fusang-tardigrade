"""
Fusang Streamlit Web Application
Web interface for Fusang phylogenetic inference tool.

Deploy to Streamlit Cloud: https://share.streamlit.io
"""

import streamlit as st
import os
import sys
import tempfile
import subprocess
import time
from pathlib import Path
from datetime import datetime

# ===========================================================
# Configuration
# ===========================================================

BASE_DIR = Path(__file__).parent
UPLOAD_FOLDER = BASE_DIR / "uploads"
RESULT_FOLDER = BASE_DIR / "static" / "results"

# Create directories
UPLOAD_FOLDER.mkdir(exist_ok=True)
RESULT_FOLDER.mkdir(parents=True, exist_ok=True)

# ===========================================================
# Page Config
# ===========================================================

st.set_page_config(
    page_title="Fusang: Tardigrade Edition",
    page_icon="🐻",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ===========================================================
# Custom CSS
# ===========================================================

st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 1rem;
    }
    .feature-box {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
    .stProgress > div > div > div > div {
        background-color: #1f77b4;
    }
</style>
""", unsafe_allow_html=True)

# ===========================================================
# Header
# ===========================================================

st.markdown('<h1 class="main-header">🐻 Fusang: Tardigrade Edition</h1>', unsafe_allow_html=True)
st.markdown("### Fast Alignment-Free Phylogenetic Inference using Spaced k-mers")
st.markdown("**IMMI Framework** | **No Alignment Required** | **Scales to 10,000+ Taxa**")
st.divider()

# ===========================================================
# Sidebar - Parameters
# ===========================================================

st.sidebar.header("⚙️ Inference Parameters")

mode = st.sidebar.selectbox(
    "Inference Mode",
    options=["auto", "l0", "l1", "l2", "l3"],
    index=0,
    help="auto: automatic level selection | l0-l3: force specific IMMI level"
)

distance_method = st.sidebar.selectbox(
    "Distance Method",
    options=["kmer", "fastme"],
    index=0,
    help="kmer: k-mer distance | fastme: FastME BIONJ+NNI"
)

kmer_k = st.sidebar.number_input(
    "k-mer size (k)",
    min_value=3,
    max_value=15,
    value=5,
    help="k-mer size for distance calculation"
)

kmer_gap = st.sidebar.selectbox(
    "Spaced Pattern (gap)",
    options=["none", "gap1", "gap2"],
    index=2,
    help="none: contiguous k-mer | gap1/gap2: spaced k-mer (indel-tolerant)"
)

tree_method = st.sidebar.selectbox(
    "Tree Method",
    options=["nj", "fastme"],
    index=0,
    help="NJ: Neighbor-Joining (O(n³)) | FastME: BIONJ+NNI (O(n²))"
)

max_group = st.sidebar.number_input(
    "Max Group Size",
    min_value=50,
    max_value=500,
    value=200,
    help="Maximum taxa per group in divide-and-conquer"
)

overlap = st.sidebar.slider(
    "Overlap Ratio",
    min_value=0.0,
    max_value=0.5,
    value=0.15,
    step=0.05,
    help="Overlap ratio between groups for boundary smoothing"
)

threads = st.sidebar.number_input(
    "Threads",
    min_value=1,
    max_value=8,
    value=4,
    help="Number of parallel threads"
)

simple_mode = st.sidebar.checkbox(
    "Simple Mode (skip UPGMA pre-sorting)",
    value=False,
    help="Faster but may reduce accuracy"
)

# ===========================================================
# Main Content
# ===========================================================

col1, col2 = st.columns([2, 1])

with col1:
    st.header("📁 Upload Sequence File")
    
    uploaded_file = st.file_uploader(
        "Upload FASTA file",
        type=["fasta", "fa", "fna", "ffn", "frn"],
        help="Upload a FASTA format sequence file (max 500 MB)"
    )
    
    if uploaded_file is not None:
        # Save uploaded file
        input_path = UPLOAD_FOLDER / uploaded_file.name
        with open(input_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
        
        st.success(f"✅ File uploaded: {uploaded_file.name}")
        
        # Show file info
        file_size = len(uploaded_file.getvalue()) / 1024 / 1024  # MB
        st.info(f"📊 File size: {file_size:.2f} MB")
        
        # Count sequences
        with open(input_path) as f:
            seq_count = sum(1 for line in f if line.startswith(">"))
        st.info(f"🧬 Sequences: {seq_count}")
        
        # Run button
        if st.button("🚀 Run Fusang", type="primary", use_container_width=True):
            with st.spinner("Running Fusang inference..."):
                # Build command
                output_path = RESULT_FOLDER / f"{input_path.stem}.nwk"
                cmd = [
                    sys.executable,
                    str(BASE_DIR / "fusang_v2.py"),
                    "-i", str(input_path),
                    "-o", str(output_path),
                    "-m", mode,
                    "-d", distance_method,
                    "--kmer_k", str(kmer_k),
                    "--tree_method", tree_method,
                    "--max_group", str(max_group),
                    "--overlap", str(overlap),
                    "-t", str(threads),
                ]
                if kmer_gap != "none":
                    cmd.extend(["--kmer_gap", kmer_gap])
                if simple_mode:
                    cmd.append("--simple")
                
                st.code(" ".join(cmd), language="bash")
                
                # Run subprocess
                progress_bar = st.progress(0)
                log_expander = st.expander("📋 Run Log", expanded=True)
                
                try:
                    process = subprocess.Popen(
                        cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        cwd=str(BASE_DIR)
                    )
                    
                    log_lines = []
                    progress = 0
                    
                    for line in iter(process.stdout.readline, ''):
                        if line:
                            log_lines.append(line.strip())
                            if any(kw in line for kw in ["[MAIN]", "[kmer]", "[DCM]", "[NNI]"]):
                                progress = min(95, progress + 2)
                                progress_bar.progress(progress)
                    
                    process.wait()
                    progress_bar.progress(100)
                    
                    with log_expander:
                        st.text("\n".join(log_lines[-50:]))
                    
                    if process.returncode == 0:
                        st.success("✅ Inference completed successfully!")
                        
                        # Read and display Newick tree
                        if output_path.exists():
                            with open(output_path) as f:
                                newick_str = f.read().strip()
                            
                            st.subheader("🌳 Result Tree (Newick)")
                            st.code(newick_str, language="text")
                            
                            # Download button
                            st.download_button(
                                label="📥 Download Tree File",
                                data=newick_str,
                                file_name=f"{input_path.stem}.nwk",
                                mime="text/plain",
                            )
                    else:
                        st.error(f"❌ Fusang failed with return code {process.returncode}")
                        with log_expander:
                            st.text("\n".join(log_lines))
                
                except Exception as e:
                    st.error(f"❌ Error: {str(e)}")

with col2:
    st.header("ℹ️ About Fusang")
    
    st.markdown("""
    <div class="feature-box">
    <b>Fusang: Tardigrade Edition</b> is a fast, alignment-free phylogenetic 
    inference tool using spaced k-mers and the IMMI framework.
    </div>
    """, unsafe_allow_html=True)
    
    st.subheader("🎯 Key Features")
    st.markdown("""
    - ✅ No sequence alignment required
    - ✅ Handles 10,000+ taxa (~70s, 609 MB RAM)
    - ✅ IMMI framework (L0-L3) with auto-selection
    - ✅ Multi-k ensemble (k=5,7,9)
    - ✅ Outperforms IQ-TREE2 GTR by 1.8× on indel-rich data
    """)
    
    st.subheader("📊 Performance")
    st.markdown("""
    | n (taxa) | Time | RAM |
    |-----------|------|-----|
    | 200 | 4.6s | <1 GB |
    | 1,000 | 18s | ~2 GB |
    | 5,000 | 42s | ~4 GB |
    | 10,000 | ~70s | ~609 MB |
    """)
    
    st.subheader("📖 Citation")
    st.markdown("""
    Kong L, Zhang L. Fast alignment-free phylogenetic inference 
    using spaced k-mers and the IMMI framework. 
    *Nucleic Acids Research* (under review), 2026.
    """)
    
    st.subheader("🔗 Links")
    st.markdown("""
    - [GitHub Repository](https://github.com/zhanglknt/fusang-tardigrade)
    - [Documentation](https://github.com/zhanglknt/fusang-tardigrade#readme)
    - [Issue Tracker](https://github.com/zhanglknt/fusang-tardigrade/issues)
    """)

# ===========================================================
# Footer
# ===========================================================

st.divider()
st.markdown(
    "Fusang: Tardigrade Edition | "
    "IMMI Framework | "
    "[GitHub](https://github.com/zhanglknt/fusang-tardigrade) | "
    "DOI: [10.5281/zenodo.20746742](https://doi.org/10.5281/zenodo.20746742)"
)
