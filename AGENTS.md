# AGENTS.md

This file defines repository-specific rules for coding agents working in this project.

## Scope

- Audience: coding agents only.
- Goal: make safe, minimal, correct changes to the Fortran library and related build/test files.

## Project Facts

- Primary language: Fortran (modernized ODRPACK95 codebase).
- Main source directory: `src/`.
- Tests directory: `test/`.
- Examples directory: `example/`.
- Build systems present: Meson, fpm, CMake.
- Preferred default validation in this repository: Meson tests.

## Hard Rules

- Make the smallest possible change that satisfies the request.
- Do not refactor unrelated code.
- Do not change public behavior unless explicitly requested.
- Keep existing coding style and naming conventions.
- Do not introduce new dependencies unless necessary and justified.
- Do not modify generated site artifacts under `_site/` unless explicitly asked.
- Never run destructive git commands (`git reset --hard`, `git checkout --`, history rewrites) unless explicitly requested.
- Do not revert user changes that are unrelated to your task.
- If unexpected unrelated modifications appear while working, stop and ask how to proceed.

## Edit Workflow

1. Read the relevant source and tests first.
2. Implement minimal edits in the most relevant file(s).
3. Add or update tests when behavior changes or bugs are fixed.
4. Run validation commands (see below).
5. Report what changed, why, and any remaining risks.

## Validation Commands

Run from repository root.

- Setup (if needed):
  - `meson setup builddir -Dbuild_tests=true`
- Build:
  - `meson compile -C builddir`
- Required default test command after edits:
  - `meson test -C builddir`

Optional secondary checks when relevant:

- `fpm test --profile debug`
- `ctest --test-dir builddir` (if using CMake workflow)

## Fortran-Specific Guidance

- Preserve numerical behavior unless a change request explicitly targets algorithms.
- Be careful with array shapes, intents, and kind usage.
- Keep interfaces consistent across modules (`odrpack_core`, `odrpack`, C API bindings).
- Prefer clear, explicit code over clever transformations in numerical kernels.

## Documentation and Examples

- If public APIs or behavior change, update relevant docs/comments and examples.
- Keep README command examples valid when changing build or test workflow.

## Completion Checklist

- Changes are minimal and focused.
- Relevant tests exist and pass with `meson test -C builddir`.
- No unrelated files were modified.
- Any limitations or follow-up work are clearly reported.
