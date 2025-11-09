import streamlit as st
import os
import tempfile
import pandas as pd
from pathlib import Path
import subprocess

st.set_page_config(page_title="DeeProMic Web", page_icon="🧬", layout="wide")
st.title("🧬 DeeProMic: Therapeutic Protein Classifier")

# Upload section
uploaded_file = st.file_uploader("Upload Protein FASTA File", type=['fasta', 'fa', 'txt'])

if uploaded_file:
    # Create temp directory
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Save uploaded file
        input_fasta = temp_path / "input.fasta"
        with open(input_fasta, "wb") as f:
            f.write(uploaded_file.getvalue())
        
        # Optional threshold
        threshold = st.slider("Probability Threshold", 0.0, 1.0, 0.5, 0.05)
        
        # Run button
        if st.button("Run Analysis", type="primary"):
            with st.spinner("Running DeeProMic analysis... This may take several minutes."):
                try:
                    # Run your existing script
                    output_dir = temp_path / "results"
                    cmd = [
                        "python", "deepromic.py",
                        "-i", str(input_fasta),
                        "-t", str(threshold),
                        "-o", str(output_dir)
                    ]
                    
                    # Execute
                    result = subprocess.run(cmd, capture_output=True, text=True, cwd="/home/beast/depromic_output")
                    
                    if result.returncode != 0:
                        st.error(f"Error: {result.stderr}")
                        st.stop()
                    
                    # Show results
                    st.success("Analysis Complete!")
                    
                    # Display and provide downloads
                    tab1, tab2, tab3, tab4 = st.tabs([
                        "📊 Probability Scores", 
                        "🎯 Potential Targets", 
                        "🔍 BLAST Results", 
                        "📥 All Downloads"
                    ])
                    
                    # Tab 1: Probability scores
                    with tab1:
                        prob_file = output_dir / "probability_score.csv"
                        if prob_file.exists():
                            df = pd.read_csv(prob_file)
                            st.dataframe(df, use_container_width=True)
                            st.download_button(
                                "Download Probability Scores",
                                df.to_csv(index=False),
                                "probability_scores.csv",
                                "text/csv"
                            )
                    
                    # Tab 2: Filtered targets
                    with tab2:
                        filtered_file = output_dir / "filtered_sequences.csv"
                        fasta_file = output_dir / "potential_targets.fasta"
                        
                        if filtered_file.exists():
                            df_filtered = pd.read_csv(filtered_file)
                            st.dataframe(df_filtered, use_container_width=True)
                            st.download_button(
                                "Download Filtered Table",
                                df_filtered.to_csv(index=False),
                                "filtered_sequences.csv",
                                "text/csv"
                            )
                        
                        if fasta_file.exists():
                            with open(fasta_file) as f:
                                st.download_button(
                                    "Download Potential Targets (FASTA)",
                                    f.read(),
                                    "potential_targets.fasta",
                                    "text/fasta"
                                )
                    
                    # Tab 3: BLAST results
                    with tab3:
                        blast_essential = output_dir / "blast_against_essential_genes.tsv"
                        blast_host = output_dir / "blast_against_host.tsv"
                        
                        if blast_essential.exists():
                            st.subheader("BLAST vs Essential Genes")
                            st.dataframe(pd.read_csv(blast_essential, sep='\t'), use_container_width=True)
                        
                        if blast_host.exists():
                            st.subheader("BLAST vs Host")
                            st.dataframe(pd.read_csv(blast_host, sep='\t'), use_container_width=True)
                    
                    # Tab 4: All files zip
                    with tab4:
                        import zipfile
                        import io
                        
                        zip_buffer = io.BytesIO()
                        with zipfile.ZipFile(zip_buffer, 'w') as zf:
                            for file in output_dir.glob("*"):
                                if file.is_file():
                                    zf.write(file, file.name)
                        
                        st.download_button(
                            "📦 Download All Results (ZIP)",
                            zip_buffer.getvalue(),
                            "deepromic_results.zip",
                            "application/zip"
                        )
                        
                except Exception as e:
                    st.error(f"An error occurred: {str(e)}")
