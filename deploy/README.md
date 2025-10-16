# PAM Deployment Scripts

This directory contains scripts for building precompiled system images for PAM.

## Scripts

### `build_sysimage.jl`

Builds a precompiled Julia system image containing PAM and its dependencies for faster startup.

**Usage:**

```bash
# Build with native CPU optimization (recommended)
julia --project=. deploy/build_sysimage.jl

# Build with generic CPU target (portable across systems)
julia --project=. deploy/build_sysimage.jl --cpu-target=generic
```

**What it does:**

1. Activates PAM project environment
2. Installs PackageCompiler if needed
3. Creates precompile script that runs PAM tests
4. Builds system image (~5-10 minutes)
5. Creates launcher script for easy usage

**Output:**

- `sysimage/sys_pam.dylib` (or `.so` on Linux) - Precompiled system image
- `sysimage/precompile_script.jl` - Script used for precompilation
- `sysimage/pam` - Launcher script

## Using the Precompiled System Image

### Method 1: Direct Julia invocation

```bash
julia --project=. --sysimage=sysimage/sys_pam.dylib
```

### Method 2: Using launcher script

```bash
# Interactive REPL
./sysimage/pam -i

# Run a script
./sysimage/pam test/PAM_test.jl examples/template_D3D_1layer_2species.json
```

## Performance Benefits

- **Normal Julia startup**: ~50 seconds for first PAM use
- **With sysimage**: ~5 seconds for first PAM use
- **Package loading**: Near-instant (precompiled)

## CPU Targets

| Target | Description | Use Case |
|--------|-------------|----------|
| `native` | Optimized for current CPU | Best performance, local use only |
| `generic` | Works on all CPUs | Portable, slightly slower |
| `apple-m1` | Apple Silicon optimized | M1/M2/M3 Macs |
| `haswell` | Intel Haswell+ | Intel CPUs 2013+ |

## Troubleshooting

**Build fails with "out of memory":**
```bash
julia --project=. deploy/build_sysimage.jl --cpu-target=generic
```

**Permission denied on launcher:**
```bash
chmod +x sysimage/pam
```

**Want to rebuild:**
```bash
rm -rf sysimage/
julia --project=. deploy/build_sysimage.jl
```

## Based On

This script is inspired by FUSE's `install_pam.jl` but simplified for PAM's needs:
- No environment management (uses current project)
- No IJulia kernel installation
- No Makefile parsing
- Minimal dependencies
