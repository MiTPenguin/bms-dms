#===============================================================================
# ubuntu base image
#===============================================================================

FROM ubuntu:24.04

#===============================================================================
# *nix utilities
#===============================================================================

RUN apt-get update && \
    apt-get install -y \
        parallel \
        wget \
        curl \
        pigz \
        vim \
        less \
        tre-agrep \
        bzip2 \
        ca-certificates \
        cargo \
        curl \
        libc6 \
        libxml2 \
        libxml2-dev \
        gcc \
        g++ \
        cmake \
        git && \
    apt-get clean

#===============================================================================
# nix environment
#===============================================================================

WORKDIR /environment

COPY flake.nix flake.nix

COPY flake.lock flake.lock

RUN curl --proto '=https' \
    --tlsv1.2 \
    -sSf \
    -L https://github.com/DeterminateSystems/nix-installer/releases/download/v0.32.3/nix-installer.sh | \
        sh -s -- install linux \
        --extra-conf "sandbox = false" \
        --init none \
        --no-confirm

ENV PATH="${PATH}:/nix/var/nix/profiles/default/bin"

RUN nix develop --command echo

CMD ["/bin/bash"]
