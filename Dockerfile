# EUROSEM web service (public mode).  Build context = the folder that contains
# headless/ and eurosem_web/ (not the deploy folder).  See deploy/README.md.
#
#   docker build -t eurosem .
#   docker run --rm -p 8080:8080 eurosem        then open http://localhost:8080

# ---- stage 1: compile the model (same flags as eurosem_web/build_exe.sh) ----
FROM python:3.12-slim-bookworm AS build
RUN apt-get update && apt-get install -y --no-install-recommends gfortran \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY headless/stage4_rval/*.for ./
RUN gfortran -std=legacy -fno-automatic -fcheck=bounds -g -w -o /src/eurosem *.for

# ---- stage 2: the service -------------------------------------------------
FROM python:3.12-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends libgfortran5 \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --system --create-home --uid 10001 eurosem
WORKDIR /app
COPY --from=build /src/eurosem /app/eurosem
COPY eurosem_web/server.py eurosem_web/app.html eurosem_web/eurosem_io.py eurosem_web/eurosem_form.py \
     eurosem_web/eurosem_run.py eurosem_web/eurosem_export.py eurosem_web/eurosem_guide.py \
     eurosem_web/eurosem_sens.py /app/
# files copied from some machines are readable only by their owner; the service runs as a different user
RUN chmod -R a+rX /app && chmod a+rx /app/eurosem
ENV EUROSEM_PUBLIC=1 EUROSEM_EXE=/app/eurosem PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 PORT=8080
USER eurosem
EXPOSE 8080
HEALTHCHECK --interval=60s --timeout=5s CMD python -c "import os,urllib.request as u;u.urlopen('http://127.0.0.1:%s/healthz'%os.environ.get('PORT','8080'),timeout=4)"
CMD ["python", "server.py"]
