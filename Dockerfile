# syntax=docker/dockerfile:1.7

# Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
# Distributed under the terms of the Apache License, Version 2.0.

# ---------- Build Stage ----------
FROM debian:bookworm-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN git clone --depth 1 https://github.com/hillerlab/CESAR2.0.git \
    && cd CESAR2.0 \
    && make

# ---------- Runtime Stage ----------
FROM debian:bookworm-slim

# perl   -> needed by the helper scripts under tools/
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        perl \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/CESAR2.0 /opt/CESAR2.0

ENV PATH="/opt/CESAR2.0:/opt/CESAR2.0/tools:${PATH}"

WORKDIR /data
CMD ["bash"]
