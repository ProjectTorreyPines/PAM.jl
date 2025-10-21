# PAM Deployment

Installation scripts for PAM across different platforms.

## Quick Start

```bash
# Generic installation (macOS, Linux, HPC)
./deploy/generic/install_pam.sh

# Perlmutter deployment (NERSC)
./deploy/perlmutter/deploy.sh

# Dev: Build sysimage only
julia --project=. deploy/generic/build_sysimage.jl
```

## Directory Structure

```
deploy/
├── generic/              # Reusable installation scripts
│   ├── install_pam.sh   # Main installer
│   ├── build_sysimage.jl
│   └── README.md        # Generic installation guide
│
└── perlmutter/          # NERSC Perlmutter specific
    ├── deploy.sh        # Main deployment
    └── README.md        # Perlmutter guide
```

## Installation Methods

### Standalone Installation (`install_pam.sh`)

Creates isolated PAM environment with optimized sysimage.

```bash
# Basic (installs to ~/.local/pam)
./deploy/generic/install_pam.sh

# Custom location
./deploy/generic/install_pam.sh --install-dir ~/pam

# Portable build
./deploy/generic/install_pam.sh --cpu-target=generic
```

**Options:**
- `--install-dir DIR` - Installation directory (default: `~/.local/pam`)
- `--cpu-target TARGET` - CPU optimization (default: `generic;native`)
- `--verbose` - Show full output

**Usage after install:**
```bash
~/.local/pam/bin/pam              # REPL
~/.local/pam/bin/pam script.jl    # Run script
```

### Sysimage Build Only (`build_sysimage.jl`)

Build optimized sysimage for existing PAM environment.

```bash
julia --project=. deploy/generic/build_sysimage.jl
```

**Options:**
- `--cpu-target=TARGET` - CPU optimization (default: `generic;native`)
- `--outdir=DIR` - Output directory (default: `sysimage/`)

## CPU Targets

| Target | Use Case |
|--------|----------|
| `generic;native` | **Default** - Portable + optimized |
| `generic` | Max compatibility |
| `native` | Single machine only |

## Platform Guides

- [Generic Installation](generic/README.md) - Details for generic/HPC installations
- [Perlmutter Deployment](perlmutter/README.md) - NERSC Perlmutter specifics

## Troubleshooting

```bash
# CPU target errors
./deploy/generic/install_pam.sh --cpu-target=generic

# Rebuild from scratch
rm -rf ~/.local/pam
./deploy/generic/install_pam.sh
```
