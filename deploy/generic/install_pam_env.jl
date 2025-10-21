#!/usr/bin/env julia
# PAM Environment Installation Script
# Creates isolated Julia environment with PAM sysimage
# Called by install_pam.sh - do not run directly

# Verify required environment variables
@assert ("PAM_INSTALL_DIR" in keys(ENV)) "Error: PAM_INSTALL_DIR environment variable not set"
@assert ("PAM_REPO_DIR" in keys(ENV)) "Error: PAM_REPO_DIR environment variable not set"

install_dir = ENV["PAM_INSTALL_DIR"]
pam_repo = ENV["PAM_REPO_DIR"]
cpu_target = get(ENV, "JULIA_CPU_TARGET", "native")

println("### PAM Environment Installation")
println("    Install directory: $install_dir")
println("    PAM repository:    $pam_repo")
println("    CPU target:        $cpu_target")
println("")

import Pkg
import TOML
import Dates

# ===== Parse PAM dependencies from Project.toml =====
println("### Reading PAM dependencies from Project.toml")
function get_pam_dependencies()
    project_file = joinpath(pam_repo, "Project.toml")
    project = TOML.parsefile(project_file)

    # Get all dependencies
    deps = get(project, "deps", Dict())
    dep_names = collect(keys(deps))

    println("    PAM dependencies: ", dep_names)
    return dep_names
end

pam_deps = get_pam_dependencies()

# ===== Setup Isolated Environment =====
println("")
println("### Activating isolated environment")
Pkg.activate(install_dir)

# Add registries (General first, then FuseRegistry)
println("### Adding Julia registries")
println("    Adding General registry...")
try
    Pkg.Registry.add("General")
catch e
    println("    General registry already added")
end

println("    Adding FuseRegistry...")
try
    Pkg.Registry.add(Pkg.RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
catch e
    println("    FuseRegistry already added")
end
Pkg.Registry.update()

# ===== Install Dependencies =====
println("")
println("### Installing PAM dependencies")

# Add PackageCompiler first (needed for sysimage)
println("    Installing PackageCompiler...")
Pkg.add("PackageCompiler")

# Add PAM's dependencies explicitly to Project.toml FIRST
# (Must be before develop, otherwise they only go into Manifest.toml)
# (PackageCompiler requires packages to be in Project.toml, not just Manifest.toml)
println("    Adding PAM dependencies explicitly for sysimage...")
Pkg.add(pam_deps)

# Add PAM itself
println("    Installing PAM package...")

# Check installation source
install_source = get(ENV, "PAM_INSTALL_SOURCE", "registry")
git_ref = get(ENV, "PAM_GIT_REF", "")

if install_source == "git" && !isempty(git_ref)
    # Install from git repository at specific branch/commit
    println("    Installing PAM from git: $git_ref")
    pam_repo_url = "https://github.com/ProjectTorreyPines/PAM.jl.git"

    try
        Pkg.add(url=pam_repo_url, rev=git_ref)
        println("    ✓ PAM installed from git ref: $git_ref")
    catch e
        println("    ✗ Failed to install from git ref '$git_ref': $e")
        println("    Falling back to registry...")
        Pkg.add("PAM")
    end
else
    # Install from registry
    pam_version = get(ENV, "PAM_VERSION", "")

    if pam_version == "latest" || isempty(pam_version)
        # Install latest from registry
        println("    Installing latest PAM from registry...")
        Pkg.add("PAM")
    else
        # Try to install specific version from registry
        # Remove 'v' prefix and any suffix for version specification
        pam_version_clean = replace(pam_version, r"^v" => "", r"-.*" => "")

        try
            println("    Attempting to install PAM@$pam_version_clean from registry...")
            Pkg.add(name="PAM", version=pam_version_clean)
        catch e
            println("    ⚠️  Version $pam_version_clean not found in registry")
            println("    Installing latest PAM from registry...")
            Pkg.add("PAM")
        end
    end
end

# Instantiate to get all transitive dependencies
println("    Resolving dependencies...")
Pkg.instantiate()
Pkg.resolve()

# Precompile packages (speeds up sysimage creation)
println("")
println("### Precompiling packages")
Pkg.precompile()

# ===== Build System Image =====
println("")
println("### Building PAM system image")

using PackageCompiler

# Determine platform-specific sysimage extension
sysimage_ext = Sys.isapple() ? "dylib" : (Sys.iswindows() ? "dll" : "so")
sysimage_path = joinpath(install_dir, "sys_pam.$sysimage_ext")

# Use precompile script from PAM repository
precompile_script = joinpath(pam_repo, "deploy", "generic", "precompile_script.jl")
if !isfile(precompile_script)
    error("Precompile script not found: $precompile_script")
end

println("    Precompile script: $precompile_script")
println("    Sysimage path:     $sysimage_path")
println("    CPU target:        $cpu_target")
println("")
println("    This will take 5-10 minutes...")
println("")

# Create sysimage with PAM and its dependencies (read from Project.toml)
sysimage_packages = ["PAM"; pam_deps]
println("    Sysimage packages: ", sysimage_packages)
println("")

create_sysimage(
    sysimage_packages;
    sysimage_path=sysimage_path,
    precompile_execution_file=precompile_script,
    cpu_target=cpu_target,
    incremental=false,  # Full static compilation
)

# Make sysimage read-only to prevent accidental modification
chmod(sysimage_path, 0o555)

# ===== Freeze Environment =====
println("")
println("### Freezing environment (read-only)")

project_file = joinpath(install_dir, "Project.toml")
manifest_file = joinpath(install_dir, "Manifest.toml")

if isfile(project_file)
    chmod(project_file, 0o444)
end
if isfile(manifest_file)
    chmod(manifest_file, 0o444)
end

# ===== Create Environment Info =====
println("")
println("### Creating environment metadata")

info_file = joinpath(install_dir, "PAM_ENV_INFO.txt")
info_content = """
PAM Environment Information
===========================

Created: $(Dates.now())
Julia Version: $(VERSION)
PAM Repository: $pam_repo
CPU Target: $cpu_target

Installed Packages:
-------------------
$(sprint(io -> Pkg.status(io=io)))

Usage:
------
Interactive REPL:
  $(install_dir)/bin/pam -i

Run script:
  $(install_dir)/bin/pam script.jl

One-liner:
  $(install_dir)/bin/pam -e 'using PAM; ...'

Direct Julia invocation:
  export JULIA_LOAD_PATH=:$install_dir
  julia --sysimage=$sysimage_path

Environment Variables:
----------------------
JULIA_DEPOT_PATH=~/.julia:$(install_dir)/.julia:
JULIA_LOAD_PATH=:$(install_dir)

Note: User depot (~/.julia) is checked first, PAM depot as fallback.
"""

write(info_file, info_content)

# ===== Installation Complete =====
println("")
println("✓ PAM environment installation complete!")
println("")
