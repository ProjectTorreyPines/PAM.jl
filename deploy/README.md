# PAM Deployment Scripts

This directory contains scripts for building precompiled system images for PAM.

## Scripts

### `build_sysimage.jl`

Builds a precompiled Julia system image containing PAM and its dependencies for faster startup.

**Usage:**

```bash
julia --project=. deploy/build_sysimage.jl [--cpu-target=TARGET] [--output=DIR]
```

**Options:**
- `--cpu-target=TARGET` - CPU optimization (default: native)
- `--output=DIR` - Output directory (default: sysimage/)
- `--help` - Show help

**Output:** `<output_dir>/sys_pam.{dylib|so}`, launcher script

## Usage

```bash
# Direct
julia --project=. --sysimage=sysimage/sys_pam.dylib

# Launcher script
./sysimage/pam -i
./sysimage/pam script.jl
```

**Performance:** ~10x faster package loading (50s → 5s)

## Troubleshooting

```bash
# Out of memory → use generic target
julia --project=. deploy/build_sysimage.jl --cpu-target=generic

# Permission denied
chmod +x sysimage/pam

# Rebuild
rm -rf sysimage/ && julia --project=. deploy/build_sysimage.jl
```
