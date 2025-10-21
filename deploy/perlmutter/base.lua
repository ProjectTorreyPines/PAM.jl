-- PAM Base Module Template for Perlmutter
-- Variables basedir, pam_version, cpu_target are set by installer

whatis("Name: PAM")
whatis("Version: " .. pam_version)
whatis("Description: Pellet Ablation Model - Julia environment")
whatis("URL: https://github.com/ProjectTorreyPines/PAM.jl")

-- Conflict with other PAM versions
family("pam")

-- Environment setup
local env_dir = pathJoin(basedir, "environments", pam_version)
local base_depot = pathJoin(env_dir, ".julia")

setenv("PAM_HOME", basedir)
setenv("PAM_VERSION", pam_version)
setenv("PAM_ENV_DIR", env_dir)

-- Julia configuration
setenv("JULIA_CPU_TARGET", cpu_target)

-- Julia depot path: Always use user's HOME depot for isolation
-- Each user gets independent package management in $HOME/.julia
-- PAM depot is appended as fallback for sysimage metadata
local home = os.getenv("HOME")
local user_depot = pathJoin(home, ".julia")
setenv("JULIA_DEPOT_PATH", user_depot .. ":" .. base_depot .. ":")

-- PAM sysimage environment is the last place Julia looks for packages
setenv("JULIA_LOAD_PATH", ":" .. env_dir)

-- Add PAM executable to PATH
prepend_path("PATH", pathJoin(env_dir, "bin"))

-- Display load message
if (mode() == "load") then
    LmodMessage("PAM " .. pam_version .. " environment loaded")
    LmodMessage("Run 'pam' to start interactive REPL")
end
