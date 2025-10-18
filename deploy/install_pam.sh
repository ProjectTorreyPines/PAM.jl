#!/bin/bash
# PAM Standalone Installation Script
# Creates isolated Julia environment with PAM sysimage
# Usage: ./deploy/install_pam.sh [OPTIONS]

set -e  # Exit on error

# ===== Configuration =====
DEFAULT_INSTALL_DIR="$HOME/.local/pam"
DEFAULT_CPU_TARGET="generic"  # Portable by default, use --cpu-target=native for max performance
PAM_REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ===== Argument Parsing =====
INSTALL_DIR="$DEFAULT_INSTALL_DIR"
CPU_TARGET="$DEFAULT_CPU_TARGET"
JULIA_THREADS=1
VERBOSE=false

usage() {
    cat << EOF
PAM Installation Script

Usage: $0 [OPTIONS]

Options:
  --install-dir DIR     Installation directory (default: $DEFAULT_INSTALL_DIR)
  --cpu-target TARGET   CPU optimization target (default: generic)
                        - generic: Portable, works on all CPUs (default)
                        - native: Fastest on your CPU, not portable
                        - haswell: Intel Haswell+ optimized
                        - znver3: AMD EPYC optimized
  --threads N           Julia threads for installation (default: 1)
  --verbose            Enable verbose output
  --help               Show this help message

Examples:
  # Basic installation (portable, works everywhere)
  $0

  # Maximum performance (only for your CPU)
  $0 --cpu-target=native

  # Custom installation directory
  $0 --install-dir ~/pam_env

  # HPC cluster optimized build
  $0 --cpu-target="generic;znver3" --install-dir /scratch/pam

EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --install-dir)
            INSTALL_DIR="$2"
            shift 2
            ;;
        --cpu-target)
            CPU_TARGET="$2"
            shift 2
            ;;
        --threads)
            JULIA_THREADS="$2"
            shift 2
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --help|-h)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# ===== Environment Detection =====
echo "===================================="
echo "PAM Installation"
echo "===================================="
echo ""

# Detect Julia
if command -v julia >/dev/null 2>&1; then
    JULIA_BIN=$(command -v julia)
    JULIA_VERSION=$(julia -e 'print(VERSION)')
    echo "✓ Found Julia $JULIA_VERSION at $JULIA_BIN"
else
    echo "✗ Julia not found in PATH"
    echo ""
    echo "Please install Julia first:"
    echo "  - Via juliaup: curl -fsSL https://install.julialang.org | sh"
    echo "  - Or download from: https://julialang.org/downloads/"
    exit 1
fi

# Detect OS
OS_TYPE=$(uname -s)
case "$OS_TYPE" in
    Linux*)
        OS_NAME="Linux"
        SYSIMAGE_EXT="so"
        ;;
    Darwin*)
        OS_NAME="macOS"
        SYSIMAGE_EXT="dylib"
        ;;
    *)
        echo "✗ Unsupported OS: $OS_TYPE"
        exit 1
        ;;
esac
echo "✓ Detected OS: $OS_NAME"

# Resolve installation directory
INSTALL_DIR=$(cd "$(dirname "$INSTALL_DIR")" 2>/dev/null && pwd)/$(basename "$INSTALL_DIR") || INSTALL_DIR=$(realpath "$INSTALL_DIR" 2>/dev/null || echo "$INSTALL_DIR")

echo ""
echo "Configuration:"
echo "  Install directory: $INSTALL_DIR"
echo "  CPU target:        $CPU_TARGET"
echo "  Julia threads:     $JULIA_THREADS"
echo "  Sysimage ext:      $SYSIMAGE_EXT"
echo ""

# ===== Pre-installation Checks =====
if [[ -d "$INSTALL_DIR" ]]; then
    read -p "Directory $INSTALL_DIR already exists. Remove and reinstall? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing installation..."
        rm -rf "$INSTALL_DIR"
    else
        echo "Installation cancelled."
        exit 0
    fi
fi

# ===== Create Installation Directory Structure =====
echo "Creating installation directory..."
mkdir -p "$INSTALL_DIR"
mkdir -p "$INSTALL_DIR/.julia"  # Isolated Julia depot
mkdir -p "$INSTALL_DIR/bin"

# ===== Set Environment Variables for Installation =====
export PAM_INSTALL_DIR="$INSTALL_DIR"
export PAM_REPO_DIR="$PAM_REPO_DIR"
export JULIA_DEPOT_PATH="$INSTALL_DIR/.julia:"
export JULIA_CPU_TARGET="$CPU_TARGET"
export JULIA_NUM_THREADS="$JULIA_THREADS"
export JULIA_PKG_USE_CLI_GIT=true  # Better for HPC environments

# ===== Run Julia Installation Script =====
echo ""
echo "===================================="
echo "Running Julia installation..."
echo "===================================="
echo ""

if [[ "$VERBOSE" == true ]]; then
    julia --startup-file=no "$PAM_REPO_DIR/deploy/install_pam_env.jl"
else
    # Show progress: major steps, package operations, and errors
    julia --startup-file=no "$PAM_REPO_DIR/deploy/install_pam_env.jl" 2>&1 | grep -E "^(###|✓|✗|Error|Warning|    |Precompiling|Downloading|Installed|Updating|Cloning|Building)" || true
fi

JULIA_EXIT_CODE=${PIPESTATUS[0]}

if [[ $JULIA_EXIT_CODE -ne 0 ]]; then
    echo ""
    echo "✗ Installation failed with exit code: $JULIA_EXIT_CODE"
    echo ""
    echo "To see full output, run with --verbose flag"
    exit $JULIA_EXIT_CODE
fi

# ===== Create User-Friendly Launcher =====
echo ""
echo "Creating launcher scripts..."

LAUNCHER_SCRIPT="$INSTALL_DIR/bin/pam"
cat > "$LAUNCHER_SCRIPT" << 'LAUNCHER_EOF'
#!/bin/bash
# PAM Launcher Script
# Auto-generated by install_pam.sh

INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SYSIMAGE_EXT="SYSIMAGE_EXT_PLACEHOLDER"

# Use user depot first (for new package installs), then PAM depot
# This follows FUSE.jl pattern: user packages go to ~/.julia, PAM packages stay isolated
export JULIA_DEPOT_PATH="$HOME/.julia:$INSTALL_DIR/.julia:"
export JULIA_LOAD_PATH=":$INSTALL_DIR"

exec julia \
    --project="$INSTALL_DIR" \
    --sysimage="$INSTALL_DIR/sys_pam.$SYSIMAGE_EXT" \
    --startup-file=no \
    "$@"
LAUNCHER_EOF

# Replace placeholders
sed -i.bak "s|SYSIMAGE_EXT_PLACEHOLDER|$SYSIMAGE_EXT|g" "$LAUNCHER_SCRIPT"
rm -f "$LAUNCHER_SCRIPT.bak"

chmod +x "$LAUNCHER_SCRIPT"

# ===== Installation Complete =====
echo ""
echo "===================================="
echo "✓ PAM Installation Complete!"
echo "===================================="
echo ""
echo "Installation location: $INSTALL_DIR"
echo ""
echo "Usage:"
echo "  1. Interactive REPL:"
echo "     $INSTALL_DIR/bin/pam"
echo ""
echo "  2. Run a script:"
echo "     $INSTALL_DIR/bin/pam myscript.jl"
echo ""
echo "  3. Verify installation:"
echo "     $INSTALL_DIR/bin/pam -e 'using PAM; println(\"PAM loaded!\")'"
echo ""
echo ""
