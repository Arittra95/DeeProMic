# Use Miniconda base image
FROM continuumio/miniconda3:latest

# Set working directory in container
WORKDIR /app

# Install system dependencies (if needed for profeatx or other tools)
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy environment file and create conda environment
COPY environment.yml .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "deepromic", "/bin/bash", "-c"]

# Install Streamlit in the environment
RUN conda install -c conda-forge streamlit

# Copy all application files to container
COPY . .

# Make profeatx executable if it's a binary/script
RUN chmod +x profeatx

# Expose Streamlit default port
EXPOSE 8501

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Run the Streamlit app
CMD ["conda", "run", "--no-capture-output", "-n", "deepromic", "streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
