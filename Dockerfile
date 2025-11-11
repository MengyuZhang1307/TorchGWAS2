# Use Ubuntu 24.04 as base
FROM ubuntu:24.04 AS builder

# Required for noninteractive installation
ENV DEBIAN_FRONTEND=noninteractive

# Install compilers and basic dependencies

RUN apt update && \
    apt install -y --no-install-recommends \
        gpg wget git gcc g++ cmake make \
        zlib1g-dev libzstd-dev dstat atop \
        python3 python3-pip python3-dev

# Install Intel MKL via official oneAPI APT repo
RUN apt-get update && apt-get install -y gnupg ca-certificates && \
    wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | \
    gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null && \
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | \
    tee /etc/apt/sources.list.d/oneAPI.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends intel-oneapi-mkl intel-oneapi-mkl-devel


# Set Intel MKL environment variables
ENV MKLROOT=/opt/intel/oneapi/mkl/latest
ENV LD_LIBRARY_PATH=${MKLROOT}/lib/intel64:${LD_LIBRARY_PATH}
ENV LIBRARY_PATH=${MKLROOT}/lib/intel64:${LIBRARY_PATH}
ENV PKG_CONFIG_PATH=${MKLROOT}/bin/mkl_link_tool

# RUN apt-get update && apt-get install -y libboost-*-dev
RUN apt-get update && apt-get install -y \
    libboost-thread-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-program-options-dev 

# Install Eigen
RUN cd /tmp && \
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz && \
tar -xf eigen-3.4.0.tar.gz && \
cp -r eigen-3.4.0/Eigen /usr/local/include/ && \
rm -rf *

# Install Armadillo
RUN cd /tmp && \
    wget https://gitlab.com/conradsnicta/armadillo-code/-/archive/14.0.1/armadillo-code-14.0.1.tar.gz && \
    tar -xf armadillo-code-14.0.1.tar.gz && \
    cp -r armadillo-code-14.0.1/include /usr/local/include/armadillo && \
    rm -rf /tmp/*
# Install SuiteSparse dependencies
RUN apt-get update && apt-get install -y \
libgmp-dev \
libmpfr-dev \
pkg-config

# Clone SuiteSparse
RUN git clone https://github.com/DrTimothyAldenDavis/SuiteSparse.git
# Create build directory
WORKDIR /SuiteSparse/build
# Configure SuiteSparse with static linking - build only essential libraries for GWAS
# Excluded: LAGraph (graph algorithms), SPEX (exact arithmetic - needs GMP), ParU (parallel LU)
# -fPIC is required for linking static libs into shared libraries (Python module)
RUN cmake .. \
    -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_STATIC_LIBS=ON \
    -DBUILD_PARU=OFF \
    -DCUDA=OFF \
    -DINSTALL_LAGRAPH_DEMOS=OFF \
    -DINSTALL_GRAPHBLAS_DEMOS=OFF \
    -DGRAPHBLAS_BUILD_TESTS=OFF \
    -DOPENMP=OFF \
    -DCMAKE_C_FLAGS="-fno-openmp -fPIC" \
    -DCMAKE_CXX_FLAGS="-fno-openmp -fPIC" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DSUITESPARSE_ENABLE_PROJECTS="suitesparse_config;mongoose;amd;btf;camd;ccolamd;colamd;cholmod;cxsparse;ldl;klu;umfpack;rbio;spqr;graphblas"

# Build everything
RUN make -j$(nproc)

# Install SuiteSparse
RUN make install


# Create a working directory
WORKDIR /GEM_BUILD

# Copy all files from your local context into the container
COPY . .

RUN echo "Files in /GEM_BUILD:" && ls -l /GEM_BUILD


# Build the project
RUN mkdir build && cd build && \
cmake .. && \
make -j$(nproc)

# Stage 2: Runtime image with CUDA support
FROM nvidia/cuda:12.6.0-runtime-ubuntu24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get clean && rm -rf /var/lib/apt/lists/*
# Install runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    time \
    python3 \
    python3-pip \
    python3-venv \
    libgomp1 \
    libboost-thread1.83.0 \
    libboost-system1.83.0 \
    libboost-filesystem1.83.0 \
    libboost-program-options1.83.0 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy MKL libraries from builder (entire directory for simplicity)
COPY --from=builder /opt/intel/oneapi/mkl/latest/lib/intel64 /opt/intel/oneapi/mkl/latest/lib/intel64

# Set MKL environment for runtime
ENV MKLROOT=/opt/intel/oneapi/mkl/latest
ENV LD_LIBRARY_PATH=${MKLROOT}/lib/intel64:/usr/local/lib:${LD_LIBRARY_PATH}
ENV MKL_THREADING_LAYER=GNU
ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

# Copy Python modules and scripts (including Mygen.so)
COPY --from=builder /GEM_BUILD/pymodules /app/pymodules
COPY --from=builder /GEM_BUILD/RunTorchGWAS.py /app/RunTorchGWAS.py

# Update library cache so the system can find the libraries
RUN ldconfig

# Create isolated Python environment and install packages
# Split into separate RUN commands and clean up after each to reduce layer size
RUN python3 -m venv /opt/venv

# Install PyTorch from CUDA index
RUN /opt/venv/bin/pip install --no-cache-dir \
    torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121 \
 && find /opt/venv -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true \
 && find /opt/venv -type f -name '*.pyc' -delete \
 && find /opt/venv -type f -name '*.pyo' -delete \
 && rm -rf /opt/venv/lib/python3.12/site-packages/torch/test

# Install other Python packages from PyPI
RUN /opt/venv/bin/pip install --no-cache-dir \
    numpy pandas scipy duckdb pyarrow tqdm \
 && find /opt/venv -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true \
 && find /opt/venv -type f -name '*.pyc' -delete \
 && find /opt/venv -type f -name '*.pyo' -delete

# Set PATH and PYTHONPATH
ENV PATH="/opt/venv/bin:$PATH"
ENV PYTHONPATH=/app:${PYTHONPATH}

WORKDIR /app

# Default to running RunTorchGWAS.py
ENTRYPOINT ["python3", "/app/RunTorchGWAS.py"]
CMD ["--help"]
