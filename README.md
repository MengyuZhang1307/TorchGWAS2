# TorchGWAS with Null Model Fitting

A GPU-accelerated GWAS analysis tool with efficient null model fitting and PyTorch-based association testing for genome-wide association studies.

## Features

- **Null Model Fitting**: CPU-based mixed model fitting with Intel MKL optimization
- **TorchGWAS**: GPU/CPU-accelerated GWAS analysis using PyTorch
- **BGEN format support**: Efficient genotype data streaming with PLINK 2.0 backend
- **Docker support**: Containerized environment with CUDA 12.6

## Requirements

### Host System Requirements
- Docker Engine (20.10+)
- NVIDIA Docker runtime (for GPU support, optional)
- NVIDIA GPU with CUDA 12.6 support (optional, can run on CPU)

**All other dependencies (GCC, CMake, MKL, Boost, etc.) are included in the Docker image!**

## Installation

### Build the Docker Image

```bash
# Clone the repository
git clone https://github.com/MengyuZhang1307/Bgen-reader-torchgwas.git
cd Bgen-reader-torchgwas

# Build the Docker image (this will take 10-15 minutes)
docker build -t torchgwas:latest .

# Verify the image was built successfully
docker images | grep torchgwas
```

**Image Details:**
- Size: ~10.8 GB
- Base: nvidia/cuda:12.6.0-runtime-ubuntu24.04
- Includes: All dependencies (GCC, CMake, MKL, Boost, PyTorch, etc.)

## Usage

### Run with GPU Support (Requires NVIDIA Docker Runtime???????)

**Note:** GPU support requires `nvidia-docker2` package installed on your host system. If not available, use CPU mode below.

```bash
docker run --gpus all \
  -v /path/to/your/data:/data \
  torchgwas:latest \
  --pheno-file /data/pheno.txt \
  --cov-file /data/cov.txt \
  --bgen /data/genotypes.bgen \
  --sample /data/samples.sample \
  --sampleid-name IID \
  --covar-names PC1 PC2 PC3 \
  --threads 10 \
  --stream-snps 1000 \
  --out /data/results.txt \
  --device cuda \
  --verbose \
  --convert
```

### Run with CPU Only (No GPU Required)

```bash
docker run --rm \
  -v /path/to/your/data:/data \
  -w /data \
  torchgwas:latest \
  --pheno-file /data/pheno.txt \
  --cov-file /data/cov.txt \
  --bgen /data/genotypes.bgen \
  --sample /data/samples.sample \
  --sampleid-name IID \
  --covar-names PC1 PC2 PC3 \
  --threads 10 \
  --stream-snps 1000 \
  --out /data/results.txt \
  --device cpu \
  --verbose \
  --convert
```

### Run Example Data

Test the installation with provided example data:

```bash
# Navigate to repository directory
cd Bgen-reader-torchgwas

# Run example analysis (CPU mode)
docker run --rm \
  -v $(pwd):/workspace \
  -w /workspace \
  torchgwas:latest \
  --pheno-file example/pheno.txt \
  --cov-file example/cov.txt \
  --bgen example/SA.bgen \
  --sample example/SA.sample \
  --sampleid-name id \
  --covar-names x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 \
  --kin-file example/kinfile.txt \
  --kin-delim tab \
  --kin-diag 1 \
  --threads 10 \
  --stream-snps 100 \
  --out example/results.txt \
  --device cpu \
  --verbose

# Results will be in example/results.txt and example/TGWAS_results.txt.parquet
```

### Interactive Shell Access

To explore the container or debug:

```bash
# With GPU
docker run --gpus all -it -v /path/to/your/data:/data torchgwas:latest /bin/bash

# CPU only
docker run -it -v /path/to/your/data:/data torchgwas:latest /bin/bash
```

## Command-Line Arguments

### Required Arguments:
- `--pheno-file`: Path to phenotype file (must have FID, IID, and phenotype columns)
- `--cov-file`: Path to covariate file (must have FID, IID, and covariate columns)
- `--bgen`: Path to BGEN genotype file
- `--sample`: Path to BGEN sample file
- `--sampleid-name`: Sample ID column header name (the **second column** in phenotype and covariate files, after FID). Examples: `IID`, `id`, `sampleid`
- `--covar-names`: Space-separated list of covariate column names to include in the model (e.g., `PC1 PC2 PC3`)

### Optional Arguments:

#### Kinship and Relatedness
- `--kin-file`: Path to kinship matrix file in three-column pairwise format (ID1, ID2, value) (default: none, assumes unrelated samples)
- `--kin-diag`: Diagonal value for kinship matrix when not accounting for inbreeding (e.g. `0.5`, default: `1.0`)

#### File Delimiters
- `--pheno-delim`: Phenotype file delimiter (default: `,`)
  - Comma: `,`
  - Tab: `tab`, `\t`, `t`, or `TAB`
  - Space: `space`, `\0`, `0`, or ` `
- `--cov-delim`: Covariate file delimiter (default: `,`)
- `--kin-delim`: Kinship file delimiter (default: `,`)

#### Variant Filtering
- `--include-snp-file`: Path to file containing subset of variants to analyze (default: all variants)
  - First line must be header: `snpid` (for PLINK or BGEN) or `rsid` (for BGEN only)
  - One variant identifier per line after header

#### Model Specification
- `--random-slope-name`: Column name in covariate file for random slope effects (default: `""`, no random slope)
- `--missing-value`: Indicator for missing values in phenotype and covariate files (default: `NA`)

#### Performance Settings
- `--threads`: Number of CPU threads for null model fitting (default: system dependent)
- `--stream-snps`: Number of SNPs to process per chunk during GWAS (default: `1000`)

#### Output Settings
- `--out`: Output file path and name (default: `out.txt`)
- `--device`: Computation device for GWAS testing (choices: `cpu` or `cuda`; default: `cuda`)
- `--verbose`: Print null model information to console (default: `False`)
- `--convert`: Convert parquet output to text file (default: `True`)
  - When `True`: Creates both `.parquet` and `.txt` output files
  - When `False`: Only creates `.parquet` output file

## Input File Formats

### Phenotype File
Tab, comma, or space-separated file with header. **Must have at least 3 columns**: FID (family ID), IID (individual ID), and at least one phenotype column. Can contain multiple phenotypes:
```
fid    id    pheno1    pheno2
1      sample1    0.5    1.2
1      sample2    0.8    1.5
2      sample3    NA     1.3
```

**Important:** 
- First column: FID (family ID)
- Second column: Individual ID - the column name must match the value provided to `--sampleid-name` (e.g., `id`, `IID`, `sampleid`)
- Remaining columns: Phenotype values
- Missing values should be coded as specified by `--missing-value` (default: `NA`)

### Covariate File
Tab, comma, or space-separated file with header. **Must have at least 3 columns**: FID (family ID), IID (individual ID), and at least one covariate column:
```
fid    id    PC1    PC2    age    sex
1      sample1    0.1    -0.2    45    1
1      sample2    0.3    0.1    52    0
2      sample3    -0.1   0.05   38    1
```

**Important:** 
- First column: FID (family ID)
- Second column: Individual ID - the column name must match the value provided to `--sampleid-name`
- Remaining columns: Covariate values
- Column names specified in `--covar-names` must exactly match the header names in this file (case-sensitive)
- Sample IDs must be in the same order as in the phenotype file

### Kinship File (optional)
Three-column pairwise format with header. Can be tab, comma, or space-separated:
```
id1    id2    kinship
sample1    sample1    1.0
sample1    sample2    0.2
sample1    sample3    0.1
sample2    sample2    1.0
sample2    sample3    0.15
sample3    sample3    1.0
```

Or with comma delimiter:
```
id1,id2,value
sample_001,sample_001,1.06823
sample_001,sample_002,0.282074
sample_002,sample_002,1.05432
```

**Notes:**
- **Three columns maximum**: ID1, ID2, and kinship value
- Header names can be any text (e.g., id1/id2/kinship or ID1/ID2/value)
- Can include all pairwise combinations or just upper/lower triangle with diagonal
- Diagonal values typically range from 0.5 to 1.0 depending on the calculation method
- Use `--kin-diag` to specify the expected diagonal value when kinship file is not provided

### Include SNP File (optional)
Single-column file with header specifying variant identifiers. First line must be either `snpid` or `rsid`:

**Using SNP IDs:**
```
snpid
chr1:100000:A:G
chr1:200000:C:T
chr2:150000:G:A
```

**Using RS IDs:**
```
rsid
rs12345
rs67890
rs11111
```

**Notes:**
- Header must be exactly `snpid` or `rsid` (case-insensitive)
- One variant identifier per line
- Used to analyze a subset of variants from the BGEN file
- All listed variants must exist in the BGEN file
- Duplicates are not allowed

## Analysis Pipeline

The analysis automatically runs two sequential stages:

### Stage 1: Null Model Fitting (CPU)
**Performed by GEMRunner using C++ backend**
- Fits a linear mixed model without genetic effects
- Accounts for population structure and relatedness using kinship matrix (if provided)
- Handles missing phenotype and covariate data 
- Uses Intel MKL with sequential threading (no OpenMP conflicts)
- Employs SuiteSparse for efficient sparse matrix operations
- Optimized with `-march=native` for your CPU architecture
- Outputs fitted null model residuals and correction factors for use in Stage 2

### Stage 2: Association Testing (GPU/CPU)
**Performed by run_gwas using PyTorch backend**
- Streams genotype data from BGEN file in chunks (`--stream-snps` parameter)
- Tests each variant for association with phenotypes using fitted null model
- Computes test statistics efficiently on GPU (CUDA) or CPU
- Supports multiple phenotypes simultaneously
- Writes results incrementally to output file

## Output

The analysis produces GWAS results in two formats:

### Parquet Output (Primary)
- **File**: `TGWAS_<outfile>.parquet`
- Binary columnar format for efficient storage and processing
- Created automatically during analysis

### Text Output (Optional)
- **File**: `<outfile>.txt` (tab-separated)
- Human-readable format
- Created by default (controlled by `--convert` flag)

### Result Columns:
- **SNP information**: 
  - rsid (variant identifier)
  - chromosome
  - position
  - alleles (reference and alternate)
- **Association statistics** (for each phenotype):
  - Beta coefficients (effect sizes)
  - Standard errors (SE)
  - T-statistics
  - P-values

**Example output format:**
```
rsid    chr    pos    ref    alt    pheno1_BETA    pheno1_SE    pheno1_T    pheno1_P    pheno2_BETA    pheno2_SE    pheno2_T    pheno2_P
rs12345    1    10000    A    G    0.05    0.02    2.5    0.012    0.03    0.015    2.0    0.045
rs67890    1    20000    C    T    -0.02    0.018    -1.1    0.27    0.01    0.012    0.83    0.40
```

## Docker Image Details

### Build Stage (ubuntu:24.04)
- **Compiler**: GCC 13 with C++20 support
- **Linear Algebra**: Intel MKL + Eigen 3.4.0 + Armadillo 14.0.1
- **Sparse Matrices**: SuiteSparse v7.8.2 (CHOLMOD, UMFPACK, SPQR)
- **I/O Libraries**: Boost (program_options, thread, system, filesystem)
- **Compression**: zstd 1.5.5, libdeflate 1.18
- **Genotype Reading**: PLINK 2.0 BGEN reader
- **Python Bindings**: pybind11 v2.12.0, fmt 11.0.2

### Runtime Stage (nvidia/cuda:12.6.0-runtime-ubuntu24.04)
- **CUDA**: 12.6 runtime for GPU acceleration (backward compatible with PyTorch CUDA 12.1)
- **PyTorch**: 2.5.1+cu121 with CUDA 12.1 support
- **MKL Runtime**: Intel MKL libraries (sequential threading, no OpenMP conflicts)
- **Python Stack**: NumPy, Pandas, SciPy, PyArrow, DuckDB, tqdm
- **Image Size**: ~10.8 GB

### Container File Locations
```
/app/                               # Application directory
├── pymodules/                      # Python modules
│   ├── Mygen.so                   # C++ extension module (GWAS core)
│   ├── ConfigueOpt.py             # Configuration wrapper
│   ├── TorchGWAS.py               # GWAS analysis functions
│   └── ReadParquet.py             # Result conversion utilities
└── RunTorchGWAS.py                # Main CLI entry point

/opt/intel/oneapi/mkl/latest/       # Intel MKL libraries
/opt/venv/                          # Python virtual environment

/workspace/                         # Default working directory (mount your data here)
```

### Key Dependencies
All libraries are statically linked into `Mygen.so` to avoid runtime dependency issues:
- fmt (v11.0.2) - String formatting
- Armadillo (v14.0.1) - Matrix operations  
- Eigen (v3.4.0) - Linear algebra
- pybind11 (v2.12.0) - Python bindings
- SuiteSparse (v7.8.2) - Sparse matrix operations

External runtime dependencies (dynamically linked):
- MKL (Intel Math Kernel Library) - BLAS/LAPACK
- Boost (v1.83.0) - thread, system, filesystem, program_options

## Performance Tips

1. **CPU Threads**: Set `--threads` to match your CPU core count for optimal null model fitting
2. **GPU Memory**: Adjust `--stream-snps` based on available VRAM:
   - 8GB GPU: `--stream-snps 500`
   - 16GB GPU: `--stream-snps 1000`
   - 24GB+ GPU: `--stream-snps 2000`
3. **Kinship Matrix**: Pre-compute and provide kinship matrix to speed up null model fitting
4. **Chunk Size**: Larger `--stream-snps` increases memory but may improve throughput
5. **Data Volume**: Always mount your data directory with `-v` for persistent results

## Troubleshooting

### Missing Required Arguments
```bash
# Error: "--covar-names is required"
# Solution: Always provide covariate names
docker run --gpus all -v /data:/data torchgwas:latest \
  --pheno-file /data/pheno.txt \
  --cov-file /data/cov.txt \
  --bgen /data/genotypes.bgen \
  --sample /data/samples.sample \
  --sampleid-name IID \
  --covar-names PC1 PC2 PC3 age sex  # <-- Required!
  
# Error: "--sampleid-name is required"
# Solution: Specify the sample ID column header name from your files
--sampleid-name id  # or IID, sampleid, etc.
```

### Column Name Mismatch
```bash
# Error: Covariate column not found
# Cause: Column names in --covar-names don't match covariate file headers
# Solution: Check your covariate file headers exactly
head -1 /data/cov.txt  # Shows: id,PC1,PC2,PC3
# Then use exact names:
--covar-names PC1 PC2 PC3  # Not pc1, not Pc1
```

### Delimiter Issues
```bash
# If you see parsing errors, verify your delimiter settings
# For tab-separated files:
--pheno-delim tab --cov-delim tab

# For space-separated files:
--cov-delim space

# For mixed delimiters:
--pheno-delim "," --cov-delim tab --kin-delim space

# Acceptable delimiter formats:
# Tab: tab, \t, t, TAB
# Space: space, \0, 0, or literal space
# Comma: ,
```

### Docker GPU Access
```bash
# Test if GPU is accessible from Docker
docker run --rm --gpus all nvidia/cuda:12.6.0-base-ubuntu24.04 nvidia-smi

# If you get error "could not select device driver with capabilities: [[gpu]]"
# You need to install nvidia-docker2:

# For Ubuntu/Debian:
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | \
  sudo tee /etc/apt/sources.list.d/nvidia-docker.list

sudo apt-get update && sudo apt-get install -y nvidia-docker2
sudo systemctl restart docker

# Verify installation
docker run --rm --gpus all nvidia/cuda:12.6.0-base-ubuntu24.04 nvidia-smi

# If GPU is not available or not needed, simply use CPU mode:
# Remove --gpus all flag and set --device cpu
```

### GPU Out of Memory
```bash
# Reduce chunk size
--stream-snps 500

# Or use CPU
--device cpu
```

### Permission Issues with Output Files
```bash
# Run with same user ID as host
docker run --gpus all --user $(id -u):$(id -g) \
  -v /path/to/your/data:/data \
  torchgwas:latest [arguments]
```

### Cannot Access Data Files
```bash
# Ensure correct path mapping
# Host path: /home/user/gwas_data/pheno.txt
# Container path: /data/pheno.txt
docker run --gpus all \
  -v /home/user/gwas_data:/data \
  torchgwas:latest \
  --pheno-file /data/pheno.txt
```

### Container Build Fails
```bash
# Clean build with no cache
docker build --no-cache -t torchgwas:latest .

# Check available disk space
df -h
```

## Example Workflow

```bash
# 1. Prepare your data directory
mkdir -p ~/gwas_analysis
cd ~/gwas_analysis

# 2. Place your input files
# - pheno.txt
# - cov.txt
# - genotypes.bgen
# - samples.sample

# 3. Build the image (one time)
docker build -t torchgwas:latest /path/to/Dockerfile/directory

# 4. Run analysis
docker run --gpus all \
  -v ~/gwas_analysis:/data \
  torchgwas:latest \
  --pheno-file /data/pheno.txt \
  --cov-file /data/cov.txt \
  --bgen /data/genotypes.bgen \
  --sample /data/samples.sample \
  --sampleid-name IID \
  --covar-names PC1 PC2 PC3 age sex \
  --threads 16 \
  --stream-snps 1000 \
  --out /data/gwas_results.txt \
  --device cuda

# 5. Results will be in ~/gwas_analysis/gwas_results.txt
```


## License

[Specify your license here]

## Contact

[Your contact information]