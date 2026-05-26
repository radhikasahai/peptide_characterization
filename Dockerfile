FROM python:3.10-slim-bookworm

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 \
    libxext6 \
    libsm6 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements-api.txt .
RUN pip install --no-cache-dir -r requirements-api.txt

COPY chem_utils.py benchmark_data.py peptide_parser.py peptide_synthesis.py peptide_visual.py ./
COPY data/benchmark_sequences.csv data/benchmark_sequences.csv
COPY api/ api/

ENV PORT=8000
EXPOSE 8000

CMD ["sh", "-c", "uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8000}"]
