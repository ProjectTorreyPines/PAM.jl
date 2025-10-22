-- PAM Base Module Template for Perlmutter
-- Placeholders will be replaced by deploy.sh

whatis("Name: PAM")
whatis("Version: _PAM_VERSION_")
whatis("Description: Pellet Ablation Model - Julia environment")
whatis("URL: https://github.com/ProjectTorreyPines/PAM.jl")

-- Load Julia module (version used during build)
depends_on("_JULIA_MODULE_")

-- Conflict with other PAM versions
family("pam")

-- Environment setup
local env_dir = "_PAM_ENV_DIR_"
local base_depot = pathJoin(env_dir, ".julia")

setenv("PAM_VERSION", "_PAM_VERSION_")
setenv("PAM_ENV_DIR", env_dir)

-- Julia configuration
setenv("JULIA_CPU_TARGET", "_CPU_TARGET_")

-- Julia depot path: user depot first, then PAM depot
-- This allows users to install their own packages while PAM environment provides base packages
local home = os.getenv("HOME")
local user_depot = os.getenv("JULIA_DEPOT_PATH")

if not user_depot then
    -- JULIA_DEPOT_PATH is unset, use default HOME user depot
    setenv("JULIA_DEPOT_PATH", pathJoin(home, ".julia") .. ":" .. base_depot .. ":")
elseif user_depot:sub(-1) == ":" then
    -- JULIA_DEPOT_PATH ends with ":", append base_depot after it
    setenv("JULIA_DEPOT_PATH", user_depot .. base_depot .. ":")
else
    -- JULIA_DEPOT_PATH doesn't end with ":", add it automatically
    setenv("JULIA_DEPOT_PATH", user_depot .. ":" .. base_depot .. ":")
end

-- PAM sysimage environment is the last place Julia looks for packages
setenv("JULIA_LOAD_PATH", ":" .. env_dir)

-- Add PAM executable to PATH
prepend_path("PATH", pathJoin(env_dir, "bin"))

-- Display load message
if (mode() == "load") then
    LmodMessage("PAM _PAM_VERSION_ environment loaded")
    LmodMessage("Run 'pam' to start interactive REPL")
end
