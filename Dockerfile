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

# Install PyTorch and Python dependencies for run_gwas
RUN pip3 install --no-cache-dir \
    torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121 \
    numpy pandas scipy

# Stage 2: Runtime image with CUDA support
FROM nvidia/cuda:12.6.0-runtime-ubuntu24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get clean && rm -rf /var/lib/apt/lists/*
# Install runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    time \
    python3 \
    python3-pip \
    libgomp1 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy MKL libraries from builder
COPY --from=builder /opt/intel/oneapi/mkl/latest/lib/intel64 /opt/intel/oneapi/mkl/latest/lib/intel64

# Set MKL environment for runtime
ENV MKLROOT=/opt/intel/oneapi/mkl/latest
ENV LD_LIBRARY_PATH=${MKLROOT}/lib/intel64:${LD_LIBRARY_PATH}
ENV MKL_THREADING_LAYER=GNU
ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

# Copy just the final binary from the builder
COPY --from=builder /GEM_BUILD/build/GEM_2.1.3 /usr/local/bin/GEM2

# Set default entry point
#ENTRYPOINT ["GEM2"]
# Copy Python modules and scripts
COPY --from=builder /GEM_BUILD/pymodules /app/pymodules
COPY --from=builder /GEM_BUILD/RunTorchGWAS.py /app/RunTorchGWAS.py

# Install Python dependencies
RUN pip3 install --no-cache-dir \
    torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121 \
    numpy pandas scipy

WORKDIR /app

# Set Python path
ENV PYTHONPATH=/app:${PYTHONPATH}

# Default to running RunTorchGWAS.py
ENTRYPOINT ["python3", "/app/RunTorchGWAS.py"]
CMD ["--help"]
