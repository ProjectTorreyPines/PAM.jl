#!/bin/bash
# PAM Perlmutter Deployment Script
# Builds PAM environment on NERSC Perlmutter compute nodes
# Based on FUSE.jl deployment pattern

set -e  # Exit on error

# ===== Parse Arguments =====
USE_LOGIN_NODE=false
PAM_BRANCH=""
PAM_COMMIT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --login)
            USE_LOGIN_NODE=true
            shift
            ;;
        --branch)
            PAM_BRANCH="$2"
            shift 2
            ;;
        --commit)
            PAM_COMMIT="$2"
            shift 2
            ;;
        --help|-h)
            cat << EOF
PAM Perlmutter Deployment Script

Usage: $0 [OPTIONS]

Options:
  --login              Build on login node (for testing only, not recommended)
  --branch BRANCH      Install from specific git branch (e.g., --branch develop)
  --commit SHA         Install from specific git commit (e.g., --commit a1b2c3d)
  --help, -h           Show this help message

Default behavior: Build on compute node via srun, install latest registry release

Version Selection Priority:
  1. --commit SHA      → Install from specific commit
  2. --branch NAME     → Install from branch HEAD
  3. (none)            → Install latest release from registry

Examples:
  # Production: Latest release from registry
  $0

  # Development: Specific branch
  $0 --branch develop

  # Bugfix: Specific commit
  $0 --commit a1b2c3d

  # Testing: Login node + development branch
  $0 --login --branch feature/new-model

Environment Variables:
  PAM_ACCOUNT    NERSC allocation (default: m4527)
  PAM_BASE_DIR   Installation base directory
EOF
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# ===== Configuration =====
# Customize these for your project/allocation
PAM_ACCOUNT="${PAM_ACCOUNT:-m4527}"  # NERSC allocation
PAM_BASE_DIR="${PAM_BASE_DIR:-/global/cfs/cdirs/m4527/codes/pam}"

# Load Julia module
command -v module >/dev/null 2>&1 || source /usr/share/lmod/lmod/init/bash

# Use specific Julia version for reproducibility
# Change this when you want to upgrade Julia version
JULIA_VERSION="${JULIA_VERSION:-1.11.4}"
module load "julia/$JULIA_VERSION"

# Verify Julia version
LOADED_JULIA=$(julia --version | grep -oP '\d+\.\d+\.\d+' || echo "unknown")
echo "Julia version: $LOADED_JULIA (required: $JULIA_VERSION)"

if [[ "$LOADED_JULIA" != "$JULIA_VERSION" ]]; then
    echo "⚠️  Warning: Loaded Julia version ($LOADED_JULIA) differs from expected ($JULIA_VERSION)"
fi

# Determine PAM version and installation source
PAM_REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$PAM_REPO_DIR"

# Priority: --commit > --branch > git tag > latest release
if [[ -n "$PAM_COMMIT" ]]; then
    PAM_VERSION="commit-$PAM_COMMIT"
    VERSION_SOURCE="git commit (specified)"
    INSTALL_SOURCE="git"
    PAM_GIT_REF="$PAM_COMMIT"
elif [[ -n "$PAM_BRANCH" ]]; then
    # Get branch HEAD commit
    if git rev-parse --verify "origin/$PAM_BRANCH" >/dev/null 2>&1; then
        BRANCH_COMMIT=$(git rev-parse --short "origin/$PAM_BRANCH")
        PAM_VERSION="branch-$PAM_BRANCH-$BRANCH_COMMIT"
        VERSION_SOURCE="git branch: $PAM_BRANCH"
        INSTALL_SOURCE="git"
        PAM_GIT_REF="$PAM_BRANCH"
    else
        echo "Error: Branch 'origin/$PAM_BRANCH' not found"
        exit 1
    fi
elif git describe --tags --exact-match >/dev/null 2>&1; then
    PAM_VERSION=$(git describe --tags --exact-match)
    VERSION_SOURCE="git tag (release)"
    INSTALL_SOURCE="registry"
    PAM_GIT_REF=""
else
    PAM_VERSION="latest"
    VERSION_SOURCE="registry (latest release)"
    INSTALL_SOURCE="registry"
    PAM_GIT_REF=""
fi

echo "===================================="
echo "PAM Perlmutter Deployment"
echo "===================================="
echo "Version:     $PAM_VERSION"
echo "Source:      $VERSION_SOURCE"
echo "Account:     $PAM_ACCOUNT"
echo "Base dir:    $PAM_BASE_DIR"
echo "Repository:  $PAM_REPO_DIR"
echo ""

if [[ "$INSTALL_SOURCE" == "git" ]]; then
    echo "Note: PAM will be installed from git repository."
    echo "      Using ref: $PAM_GIT_REF"
else
    echo "Note: PAM will be installed from FuseRegistry."
    if [[ "$PAM_VERSION" != "latest" ]]; then
        echo "      If version $PAM_VERSION exists in registry, it will be used."
        echo "      Otherwise, the latest registry version will be installed."
    fi
fi
echo ""

# ===== Directory Setup =====
ENV_DIR="$PAM_BASE_DIR/environments/$PAM_VERSION"
MODULE_DIR="$PAM_BASE_DIR/modules/pam"
BUILD_DIR="$SCRATCH/pam_build/$PAM_VERSION"
LOG_DIR="$SCRATCH/pam_build_logs"

# Check if already installed
if [[ -d "$ENV_DIR" ]]; then
    echo "✓ PAM $PAM_VERSION already installed in $ENV_DIR"
    echo ""
    read -p "Reinstall? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Deployment cancelled."
        exit 0
    fi
    echo "Removing existing installation..."
    rm -rf "$ENV_DIR"
fi

# Create directories
mkdir -p "$MODULE_DIR"
mkdir -p "$LOG_DIR"
[[ -d "$BUILD_DIR" ]] && rm -rf "$BUILD_DIR"

LOG_FILE="$LOG_DIR/$PAM_VERSION-$(date +%Y%m%d-%H%M%S).log"

echo "Build dir:   $BUILD_DIR"
echo "Log file:    $LOG_FILE"
echo ""

# ===== Detect Julia Module =====
# Extract currently loaded Julia module name (e.g., "julia/1.11.4")
JULIA_MODULE=$(module list 2>&1 | grep -oP 'julia/[\d.]+' | head -1)

if [[ -z "$JULIA_MODULE" ]]; then
    echo "✗ Error: No Julia module loaded"
    echo "  Please load Julia module first: module load julia"
    exit 1
fi

echo "Julia module: $JULIA_MODULE"
echo ""

# ===== Environment Variables for Installation =====
export PAM_INSTALL_DIR="$BUILD_DIR"
export PAM_REPO_DIR="$PAM_REPO_DIR"
export PAM_VERSION="$PAM_VERSION"
export PAM_BASE_DIR="$PAM_BASE_DIR"

# PAM installation source
export PAM_INSTALL_SOURCE="$INSTALL_SOURCE"  # "git" or "registry"
export PAM_GIT_REF="$PAM_GIT_REF"            # branch/commit if git install

# Julia configuration
export JULIA_DEPOT_PATH="$BUILD_DIR/.julia:"
export JULIA_PKG_USE_CLI_GIT=1
export JULIA_NUM_THREADS=1

# Perlmutter CPU optimization
# znver3 for AMD EPYC Milan (login and compute nodes)
# Remove problematic instructions: xsaveopt, rdrnd (like Julia official builds)
# Use clone_all for better multiversioning
export JULIA_CPU_TARGET="generic;znver3,-xsaveopt,-rdrnd,clone_all"

# Jupyter configuration (if using IJulia)
export IJULIA_NODEFAULTKERNEL=1
export JUPYTER_DATA_DIR="$BUILD_DIR/.jupyter"

echo "===================================="
echo "Environment Configuration"
echo "===================================="
echo "JULIA_CPU_TARGET:   $JULIA_CPU_TARGET"
echo "JULIA_DEPOT_PATH:   $JULIA_DEPOT_PATH"
echo "JULIA_NUM_THREADS:  $JULIA_NUM_THREADS"
echo ""

# ===== Run Installation on Compute Node =====
echo "===================================="
echo "Submitting build job to compute node"
echo "===================================="
echo ""
echo "This will:"
echo "  - Request 1 CPU node via srun"
echo "  - Time limit: 4 hours"
echo "  - Run PAM installation script"
echo "  - Build optimized system image"
echo ""
echo "Log output will be saved to: $LOG_FILE"
echo ""

# Clear any SLURM environment variables (important for scrontab)
unset ${!SLURM_@}

# ===== Run Installation =====
if [[ "$USE_LOGIN_NODE" == true ]]; then
    echo "⚠️  WARNING: Building on login node (testing mode)"
    echo "    This may be slower and CPU target may not match compute nodes"
    echo ""

    # Run directly on login node
    bash "$(dirname "${BASH_SOURCE[0]}")/../generic/install_pam.sh" \
        --install-dir "$BUILD_DIR" \
        --cpu-target "$JULIA_CPU_TARGET" \
        --threads "$JULIA_NUM_THREADS" \
        --verbose 2>&1 | tee "$LOG_FILE"
    JULIA_EXIT_CODE=${PIPESTATUS[0]}
else
    # Submit job to compute node (recommended)
    # Use srun for interactive build (better for debugging than sbatch)
    echo "Submitting job to compute node..."
    srun \
        --nodes=1 \
        --qos=regular \
        --time=03:00:00 \
        --constraint=cpu \
        --account="$PAM_ACCOUNT" \
        --output="$LOG_FILE" \
        bash "$(dirname "${BASH_SOURCE[0]}")/../generic/install_pam.sh" \
            --install-dir "$BUILD_DIR" \
            --cpu-target "$JULIA_CPU_TARGET" \
            --threads "$JULIA_NUM_THREADS" \
            --verbose

    JULIA_EXIT_CODE=$?
fi

# ===== Check Build Status =====
echo ""
if [[ $JULIA_EXIT_CODE -ne 0 ]]; then
    echo "✗ Build failed with exit code: $JULIA_EXIT_CODE"
    echo ""
    echo "Last 50 lines of log:"
    echo "===================="
    tail -n 50 "$LOG_FILE"
    echo ""
    echo "Full log: $LOG_FILE"
    exit $JULIA_EXIT_CODE
fi

# ===== Generate Module File =====
echo "===================================="
echo "Generating module file"
echo "===================================="
echo ""

MODULE_TEMPLATE="$(dirname "${BASH_SOURCE[0]}")/base.lua"
MODULE_OUTPUT="$BUILD_DIR/pam_module.lua"

echo "Reading template: $MODULE_TEMPLATE"
echo "Writing module file: $MODULE_OUTPUT"

# Substitute placeholders with actual values
sed -e "s|_PAM_ENV_DIR_|$ENV_DIR|g" \
    -e "s|_PAM_VERSION_|$PAM_VERSION|g" \
    -e "s|_CPU_TARGET_|$JULIA_CPU_TARGET|g" \
    -e "s|_JULIA_MODULE_|$JULIA_MODULE|g" \
    "$MODULE_TEMPLATE" > "$MODULE_OUTPUT"

# Validate: check for unreplaced placeholders
if grep -qE '_PAM_|_JULIA_|_CPU_' "$MODULE_OUTPUT"; then
    echo "✗ Error: Unreplaced placeholders detected in module file!"
    echo ""
    echo "Found:"
    grep -E '_PAM_|_JULIA_|_CPU_' "$MODULE_OUTPUT"
    echo ""
    echo "This indicates a mismatch between template and deploy.sh variables."
    exit 1
fi

echo "✓ Module file generated and validated"
echo ""

# ===== Deploy to Production =====
echo "===================================="
echo "Deploying to production"
echo "===================================="
echo ""

# Copy build directory to production location
echo "Copying $BUILD_DIR -> $ENV_DIR"

# Ensure parent directory exists
PARENT_DIR="$(dirname "$ENV_DIR")"
if [[ ! -d "$PARENT_DIR" ]]; then
    echo "Creating parent directory: $PARENT_DIR"
    mkdir -p "$PARENT_DIR"
fi

# Remove existing installation if present
if [[ -d "$ENV_DIR" ]]; then
    echo "Removing existing installation: $ENV_DIR"
    rm -rf "$ENV_DIR"
fi

/bin/cp -r "$BUILD_DIR" "$ENV_DIR"

# Copy module file to module directory
echo "Installing module file: $MODULE_DIR/$PAM_VERSION.lua"

# Ensure module directory exists
if [[ ! -d "$MODULE_DIR" ]]; then
    echo "Creating module directory: $MODULE_DIR"
    mkdir -p "$MODULE_DIR"
fi

/bin/cp "$BUILD_DIR/pam_module.lua" "$MODULE_DIR/$PAM_VERSION.lua"

# Update default symlink
if [[ -L "$MODULE_DIR/default.lua" ]]; then
    rm "$MODULE_DIR/default.lua"
fi
/bin/ln -s "$MODULE_DIR/$PAM_VERSION.lua" "$MODULE_DIR/default.lua"

# Verify installed PAM version
echo ""
echo "Verifying installed PAM version..."
ACTUAL_PAM_VERSION=$(julia --startup-file=no --project="$ENV_DIR" -e 'import Pkg; for pkg in values(Pkg.dependencies()); if pkg.name == "PAM"; println(pkg.version); break; end; end' 2>/dev/null || echo "unknown")

echo ""
echo "===================================="
echo "✓ PAM $PAM_VERSION Deployment Complete!"
echo "===================================="
echo ""
echo "Deployment version:  $PAM_VERSION ($VERSION_SOURCE)"
echo "Installed PAM:       v$ACTUAL_PAM_VERSION (from FuseRegistry)"
if [[ "$PAM_VERSION" != "v$ACTUAL_PAM_VERSION" ]] && [[ ! "$PAM_VERSION" =~ ^dev- ]]; then
    echo ""
    echo "⚠️  Note: Deployment version differs from installed version"
    echo "    This may occur if the registry doesn't have version $PAM_VERSION yet."
fi
echo ""
echo "Installation location: $ENV_DIR"
echo ""
echo "Usage:"
echo "  1. Load module:"
echo "     module use $PAM_BASE_DIR/modules"
echo "     module load pam/$PAM_VERSION"
echo ""
echo "  2. Run PAM REPL:"
echo "     pam"
echo ""
echo "  3. Run PAM script:"
echo "     pam my_script.jl"
echo ""
echo "  4. Check installation:"
echo "     module load pam"
echo "     pam -e 'using PAM; println(\"PAM loaded!\")'"
echo ""
