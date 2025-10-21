# PAM Perlmutter Deployment

Automated deployment for PAM on NERSC Perlmutter with AMD EPYC optimization.

## Quick Start

```bash
cd PAM.jl
./deploy/perlmutter/deploy.sh
```

Creates optimized installation with module file.

## Version Control

```bash
# Production: Latest registry release
./deploy.sh

# Development: Specific branch
./deploy.sh --branch develop

# Bugfix: Specific commit
./deploy.sh --commit a1b2c3d
```

Version priority: `--commit` > `--branch` > git tag > latest release

## Configuration

Edit `deploy.sh` for your allocation:

```bash
PAM_ACCOUNT="m4527"  # Your NERSC allocation
PAM_BASE_DIR="/global/cfs/cdirs/${PAM_ACCOUNT}/codes/pam"
```

## Architecture

```
deploy/
├── generic/                    # Reusable scripts
│   ├── install_pam.sh
│   └── install_pam_env.jl
└── perlmutter/                # Platform-specific
    ├── deploy.sh              # Main deployment (run this)
    ├── install_pam_perlmutter.jl
    └── base.lua               # Module template
```

**Design:** Thin Perlmutter wrapper delegates to generic scripts, following FUSE pattern.

## Deployment Process

1. Detect PAM version (git tag or `latest`)
2. Submit build job to compute node (`srun`)
3. Build optimized sysimage with AMD EPYC tuning
4. Copy to shared directory
5. Generate module file

## CPU Target

```bash
JULIA_CPU_TARGET="generic;znver3,-xsaveopt,-rdrnd,clone_all"
```

- `generic` - Fallback for any CPU
- `znver3` - AMD EPYC Milan optimization
- `-xsaveopt,-rdrnd` - Exclude problematic instructions
- `clone_all` - Enable multiversioning

**Benefits:** Portable (login + compute nodes) + optimized for Perlmutter CPUs.

## Directory Structure

```
/global/cfs/cdirs/m4527/codes/pam/
├── environments/
│   ├── v2.0.0/              # Tagged release
│   │   ├── sys_pam.so
│   │   ├── bin/pam
│   │   └── .julia/
│   └── branch-develop-abc/  # Development build
└── modules/pam/
    ├── v2.0.0.lua
    └── default.lua -> v2.0.0.lua
```

## Usage

```bash
# Load module
module use /global/cfs/cdirs/m4527/codes/pam/modules
module load pam

# Use PAM
pam                           # REPL
pam my_simulation.jl          # Run script
pam -e 'using PAM; ...'       # One-liner
```

Add to `~/.bashrc` for automatic access:
```bash
module use /global/cfs/cdirs/m4527/codes/pam/modules
```

## Julia Version

Pinned to **Julia 1.11.4** for sysimage compatibility.

Override if needed:
```bash
JULIA_VERSION=1.11.5 ./deploy.sh
```

## Troubleshooting

### CPU Instruction Errors

```bash
# Use generic target only
export JULIA_CPU_TARGET="generic"
./deploy.sh
```

### Build Timeout

```bash
# Edit deploy.sh, increase --time value as needed
srun --time=HH:MM:SS ...
```

### Version Mismatch

Expected when deployment name differs from installed package:
```
⚠️  Deployment: branch-develop-abc123
    Installed:  v2.0.2 (latest registry)
```

This is normal for branch/commit installations.

## Version Management

### Deploy New Version

```bash
git tag v2.1.0
git push --tags
./deploy.sh  # Auto-detects v2.1.0
```

### Switch Default

```bash
cd /global/cfs/cdirs/m4527/codes/pam/modules/pam
ln -sf v2.1.0.lua default.lua
```

### Development Build

```bash
./deploy.sh --branch feature-x  # Creates branch-feature-x-<hash>
module load pam/branch-feature-x-<hash>
```

## Maintenance

```bash
# Clean old versions
rm -rf environments/v1.*
rm modules/pam/v1.*.lua

# Update Julia version
# Edit deploy.sh: JULIA_VERSION="1.12.0"
./deploy.sh
```

## References

- [FUSE Deployment](https://github.com/ProjectTorreyPines/FUSE.jl/tree/master/deploy/perlmutter) - Pattern reference
- [Perlmutter Docs](https://docs.nersc.gov/systems/perlmutter/)
- [PackageCompiler](https://github.com/JuliaLang/PackageCompiler.jl)
