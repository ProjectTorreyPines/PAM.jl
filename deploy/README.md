# PAM Deployment Scripts

This directory contains scripts for building precompiled system images for PAM.

## Quick Start

### Standalone Installation (Recommended)

For a complete, isolated PAM environment:

```bash
./deploy/install_pam.sh
```

This creates a self-contained PAM installation with its own Julia depot and sysimage.

### Simple Sysimage Build

For development or if you just want a faster-loading PAM in your existing environment:

```bash
julia --project=. deploy/build_sysimage.jl
```

## Installation Methods

### Method 1: Standalone Installation (`install_pam.sh`)

**Use when:** Production, HPC, or sharing with others.

Features:
- Self-contained environment (works without dev repo)
- User packages go to `~/.julia`, PAM stays isolated
- Optimized sysimage (much faster startup)
- Portable launcher script

**Usage:**

```bash
# Basic installation to ~/.local/pam
./deploy/install_pam.sh

# Custom location
./deploy/install_pam.sh --install-dir ~/my_pam

# Portable build (works on any CPU)
./deploy/install_pam.sh --cpu-target=generic

# HPC optimized
./deploy/install_pam.sh --cpu-target="generic;znver3" --install-dir /scratch/$USER/pam
```

**Options:**
- `--install-dir DIR` - Installation directory (default: ~/.local/pam)
- `--cpu-target TARGET` - CPU optimization (default: **generic**)
  - `generic`: Works on all CPUs (portable, default)
  - `native`: Optimize for current CPU (faster, not portable)
  - `"generic;haswell"`: Generic + Intel Haswell+ optimized
  - `"generic;znver3"`: Generic + AMD EPYC optimized
- `--threads N` - Julia threads during installation (default: 1)
- `--verbose` - Show full output
- `--help` - Show help

**After installation:**
```bash
# Interactive REPL
~/.local/pam/bin/pam -i

# Run a script
~/.local/pam/bin/pam myscript.jl

# One-liner
~/.local/pam/bin/pam -e 'using PAM; println("Hello!")'

# Add to PATH (optional)
export PATH="$HOME/.local/pam/bin:$PATH"
```

### Method 2: Simple Sysimage Build (`build_sysimage.jl`)

**Use when:** You're developing PAM and want faster loading in your current environment.

Builds a precompiled Julia system image containing PAM and its dependencies for faster startup.

**Usage:**

```bash
julia --project=. deploy/build_sysimage.jl [--cpu-target=TARGET] [--outdir=DIR]
```

**Options:**
- `--cpu-target=TARGET` - CPU optimization (default: **generic**)
  - `generic`: Works everywhere (portable, default)
  - `native`: Fastest on your CPU (not portable)
- `--outdir=DIR` - Output directory (default: sysimage/)
- `--help` - Show help

**Output:** `<output_dir>/sys_pam.{dylib|so}`, launcher script

### `precompile_script.jl`

Traces PAM usage patterns for sysimage compilation. Test independently:
```bash
julia --project=. deploy/precompile_script.jl
```

## Performance

**Startup:** Significantly faster with sysimage (seconds vs minutes)
**Computation:** `native` target provides modest performance improvement over `generic`

## Comparison

| Feature | `install_pam.sh` | `build_sysimage.jl` |
|---------|------------------|---------------------|
| **Setup** | One command, fully automated | Requires PAM project |
| **Portability** | Works anywhere (no dev repo needed) | Tied to current environment |
| **User packages** | Install to `~/.julia` (shared) | Install to project |
| **Use case** | Production, HPC, sharing | Local development |
| **Time** | 10-15 min | 5-10 min |

## CPU Target Guide

| Use Case | Target | Why |
|----------|--------|-----|
| **Sharing/HPC** | `generic` (default) | Works everywhere |
| **Max performance** | `native` | Fastest on your CPU only |
| **HPC optimized** | `"generic;znver3"` (AMD) / `"generic;haswell"` (Intel) | Portable + optimized |

## Troubleshooting

```bash
# Out of memory during build
julia --project=. deploy/build_sysimage.jl --cpu-target=generic

# Permission denied on launcher
chmod +x sysimage/pam  # or ~/.local/pam/bin/pam

# Rebuild from scratch
rm -rf sysimage/ && julia --project=. deploy/build_sysimage.jl
rm -rf ~/.local/pam && ./deploy/install_pam.sh

# User packages not visible in PAM
# This is expected - PAM uses isolated environment
# User packages install to ~/.julia and are accessible
```
