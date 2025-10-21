# Generic PAM Installation

Platform-independent installation scripts for PAM (macOS, Linux, HPC).

## Quick Start

```bash
cd PAM.jl
./deploy/generic/install_pam.sh
```

Installs to `~/.local/pam` with optimized sysimage.

## Files

- `install_pam.sh` - Main installation script
- `install_pam_env.jl` - Julia environment setup
- `build_sysimage.jl` - Sysimage builder
- `precompile_script.jl` - Precompilation workload

## Installation Options

```bash
# Custom directory
./deploy/generic/install_pam.sh --install-dir ~/my_pam

# Portable build (no CPU optimization)
./deploy/generic/install_pam.sh --cpu-target=generic

# Maximum performance (non-portable)
./deploy/generic/install_pam.sh --cpu-target=native

# Fat binary (default, recommended)
./deploy/generic/install_pam.sh --cpu-target="generic;native"
```

**Available flags:**
- `--install-dir DIR` - Installation directory (default: `~/.local/pam`)
- `--cpu-target TARGET` - CPU optimization (default: `generic;native`)
- `--threads N` - Julia threads during build
- `--verbose` - Show full output
- `--help` - Show help

## CPU Targets

| Target | Compatibility | Performance |
|--------|--------------|-------------|
| `generic` | Universal | Baseline |
| `native` | Build CPU only | Maximum |
| `generic;native` | Universal + optimized | **Recommended** |

## Usage After Installation

```bash
# REPL
~/.local/pam/bin/pam

# Run script
~/.local/pam/bin/pam my_script.jl

# One-liner
~/.local/pam/bin/pam -e 'using PAM; println("Hello!")'

# Add to PATH (optional)
export PATH="$HOME/.local/pam/bin:$PATH"
```

## Building Sysimage Only

If PAM already installed and you only need sysimage:

```bash
julia --project=. deploy/generic/build_sysimage.jl
```

**Options:**
- `--cpu-target=TARGET` - CPU optimization
- `--outdir=DIR` - Output directory (default: `sysimage/`)

## Installation Directory Structure

```
~/.local/pam/
├── Project.toml
├── Manifest.toml
├── sys_pam.so           # Optimized sysimage
├── bin/pam              # Launcher
└── .julia/              # PAM depot (fallback for sysimage)
```

Launcher uses `$HOME/.julia` first, then PAM depot as fallback.

## Troubleshooting

```bash
# CPU instruction errors (AES-NI, xsaveopt)
./deploy/generic/install_pam.sh --cpu-target=generic

# Permission denied
chmod +x deploy/generic/install_pam.sh

# Verbose output for debugging
./deploy/generic/install_pam.sh --verbose
```

## Updating

```bash
# Remove old installation
rm -rf ~/.local/pam

# Reinstall
./deploy/generic/install_pam.sh
```

## Platform-Specific

For NERSC Perlmutter optimization, see [../perlmutter/README.md](../perlmutter/README.md).
