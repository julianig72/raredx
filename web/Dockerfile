# raredx web — containerized clinical VCF analysis service
FROM python:3.13-slim

WORKDIR /app

# system deps kept minimal; the pipeline talks to REST APIs over HTTPS
RUN pip install --no-cache-dir --upgrade pip

# install web deps
COPY web/requirements.txt /app/web/requirements.txt
RUN pip install --no-cache-dir -r /app/web/requirements.txt
ENV COPILOT_CLI_EXTRACT_DIR=/opt/copilot-runtime
RUN mkdir -p ${COPILOT_CLI_EXTRACT_DIR} \
    && python -m copilot download-runtime \
    && chmod -R a+rx ${COPILOT_CLI_EXTRACT_DIR}

# app code: the pipeline module + the web server/static assets
COPY raredx_pipeline.py /app/raredx_pipeline.py
COPY web /app/web

# non-root user (no PHI should be baked into the image; data lives in a mounted volume)
RUN useradd -m raredx && mkdir -p /data && chown raredx /data
USER raredx

ENV RAREDX_DATA_DIR=/data \
    RAREDX_MAX_MB=50 \
    RAREDX_MAX_VARIANTS=50000 \
    RAREDX_JOB_TIMEOUT_SECONDS=3600 \
    RAREDX_LLM_PROVIDER=auto \
    RAREDX_COPILOT_MODEL=gpt-5-mini \
    PORT=8000
# ENV RAREDX_EMAIL=your@institution.org   # optional NCBI E-utilities contact

EXPOSE 8000

# Reminder: this container has NO auth. Run it behind an authenticating reverse
# proxy (OAuth2/OIDC) + TLS before any clinical use.
CMD ["uvicorn", "web.server:app", "--host", "0.0.0.0", "--port", "8000"]
