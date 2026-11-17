# Reduction validation, 2026-10-26

These logs test the reduced active working tree, not the archived original sources.
Run commands from the owning repository. Julia checks used Julia 1.12.6 and the existing
local development manifests. Other Julia versions in the compatibility range were not tested.
The restricted environment used a writable temporary depot followed by the existing user
depot, with `--compiled-modules=existing`. Registry access was unavailable; installed
local dependencies sufficed. Documentation builds skipped deployment.

````sh
julia --project=test test/runtests.jl
julia --project=docs docs/make.jl
````
