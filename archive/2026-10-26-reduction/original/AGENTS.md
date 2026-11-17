# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## Project

`ClusterAlgebras.jl` is a Julia package for computing with cluster algebras: quivers and seed mutation, principal coefficients (c-/g-vectors, F-polynomials), root systems and finite-type classification, mutation classes and exchange graphs, green sequences, friezes, and Grassmannian cluster structures. It builds on `AbstractAlgebra.jl`.

**Orient first:** read `PLANNING/ARCHITECTURE.md` before starting any task. It has the file map, type hierarchy, key invariants, and mutation dispatch. The other `PLANNING/` files cover the roadmap, design spec, math background, and open tasks.

## Layout

- `src/ClusterAlgebras.jl` is the entry point: `include` order, exports, precompile workload. All other `src/` files are one topic each (see the file map in `PLANNING/ARCHITECTURE.md`).
- `ext/` holds package extensions (Graphs, GraphMakie, Latexify) gated on weak dependencies in `Project.toml`.
- `test/runtests.jl` includes one `test_*.jl` file per `src/` file; tests use `@testset` blocks.
- `examples/` contains Jupyter notebooks. Do not read `.ipynb` files - they are large; ask the user instead.

## Common commands

```sh
# Run a single test file (preferred - see workflow rules)
julia --project=. test/test_mutation.jl

# Run the full suite
julia --project=. -e "using Pkg; Pkg.test()"

# Start a REPL with the project environment active
julia --project=.
```

## Workflow rules

- **Plan before coding.** Discuss and agree on the design before writing Julia code. Once agreed, write the code directly, run the tests once, and fix from the output - no extended pre-implementation analysis.
- **Test scope.** Run only the test file(s) relevant to the change. Run the full suite only when shared infrastructure (`quiver.jl`, `seed.jl`, `mutation.jl`, exports) changes.
- **After code changes**, update the affected `PLANNING/` files directly - mark items done in `PLANNING/roadmap.md`, update `PLANNING/ARCHITECTURE.md` (public API, structure, test count) and any other affected planning docs. Do not keep a `RECENT_CHANGES.md` staging log.
- **Git:** never `git commit` or `git push` unless explicitly asked. Never add `Co-Authored-By` trailers to commit messages.

## Julia conventions

- Dependencies are declared in `Project.toml`; `Manifest.toml` is gitignored (correct for a library - consumers resolve their own versions).
- New functionality goes in the topical `src/` file (or a new one included from the entry point) with a matching `test_*.jl` included from `runtests.jl`.
- **References:** cite non-standard mathematics in `src/` comments by key from `PLANNING/references.md` (`# [FZ4] §6`), and add the entry there when you first lean on a source. Docstrings use a plain author-year instead of a key. `REFERENCES/` holds definitions extracted verbatim from the FWZ book chapters, for transcription checks.
- **Convention ledgers:** a sign, normalization, or transcription choice that is *ours* rather than forced goes in a ledger block at the top of the owning `src/` file - what was pinned, and which oracle pins it. `block_decomposition.jl` (the seven-block table and its two independent checks) is the model.
