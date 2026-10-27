#!/usr/bin/env julia
# Render every Pluto notebook in docs/notebooks/*.jl to a LEAN standalone
# HTML at docs/notebooks-static/<slug>.html via Snapshot.jl's Therapy emit
# (therapy=true): rendered cell outputs as plain HTML + WASM islands where a
# @bind group compiles — and NO embedded Pluto statefile.
#
# Cache is keyed on a SHA-256 fingerprint of the source notebook, the package
# source, Snapshot's locked source tree, and the export contract (NOT mtime, so
# `git checkout` doesn't trigger spurious rebuilds). The fingerprint is
# embedded in both exported forms; renderer/theme upgrades therefore
# invalidate cleanly.
#
# Workflow:
#   1. Render locally on whatever GPU you have (Metal / CUDA / ROCm) —
#      this populates docs/notebooks-static/<slug>.html with the fingerprint baked in.
#   2. `git add docs/notebooks-static/` and commit.
#   3. CI inherits the rendered HTML, sees a fingerprint match, skips re-render,
#      deploys what you committed. Re-render only happens when somebody
#      actually edits a .jl source (or the package source under src/).
#
# Force a full rebuild with: VASCULATURE_FORCE_NB_REBUILD=1 julia ...
#
# Usage:
#   julia --project=docs/build_env docs/extract_all.jl

import Pkg
let build_env = joinpath(@__DIR__, "build_env")
    if Base.active_project() != joinpath(build_env, "Project.toml")
        Pkg.activate(build_env; io = devnull)
    end
end

using Pluto
using Snapshot: export_notebook
using SHA: bytes2hex, sha256
using UUIDs: UUID

const NOTEBOOKS_DIR = joinpath(@__DIR__, "notebooks")
const STATIC_OUT = joinpath(@__DIR__, "notebooks-static")
const DATA_DIR = joinpath(NOTEBOOKS_DIR, "data")
const DATA_PROVENANCE = joinpath(NOTEBOOKS_DIR, "DATA_PROVENANCE.sha256")
const HASH_META = "vasculature-export-fingerprint"
const SNAPSHOT_UUID = UUID("1e64ef43-5f79-4f6b-8e97-159df4e27032")
# slug => (bond names,): @bind controls whose closures are outside WasmTarget's
# scope. They stay interactive in Pluto and publish as a disabled control.
const FORCE_FALLBACK_BONDS = Dict{String, Tuple{Vararg{Symbol}}}(
    # every control drives a CairoMakie figure
    "00_cardiac_pipeline" => (
        :z_xcat, :z_domain, :min_radius_shown, :z_flow,
        :healthy_state, :healthy_time, :healthy_z,
        :stenosis_state, :stenosis_time, :stenosis_z,
    ),
)
const EXPORT_CONTRACT = "vasculature-lean-v1|therapy=true|fragment=true|islands=true|verify=true|optimize=size|forced-fallback="

"""Discover every Pluto notebook under `docs/notebooks/`."""
function notebook_slugs()
    isdir(NOTEBOOKS_DIR) || return String[]
    return [
        splitext(f)[1]
            for f in sort(readdir(NOTEBOOKS_DIR))
            if endswith(f, ".jl") && !endswith(f, ".sessions.toml")
    ]
end

"""SHA-256 of the source .jl content."""
source_hash(slug::AbstractString) =
    bytes2hex(sha256(read(joinpath(NOTEBOOKS_DIR, "$(slug).jl"))))

"Hash Project.toml while ignoring the package's release-only version field."
function project_input_hash(path::AbstractString)
    content = read(path, String)
    normalized = replace(
        content,
        r"(?m)^version\s*=\s*\"[^\"]*\"\s*$" => "version = \"<release-version>\"",
    )
    return bytes2hex(sha256(normalized))
end

"Digest committed package code and the (optional) committed data provenance."
const _PACKAGE_INPUT_HASH = Ref{Union{Nothing, String}}(nothing)
function package_input_hash()
    _PACKAGE_INPUT_HASH[] === nothing || return _PACKAGE_INPUT_HASH[]::String
    root = dirname(@__DIR__)
    project = joinpath(root, "Project.toml")
    files = isfile(DATA_PROVENANCE) ? [DATA_PROVENANCE] : String[]
    for dir in (joinpath(root, "src"),)
        isdir(dir) || continue
        for (base, _, names) in walkdir(dir), name in names
            path = joinpath(base, name)
            isfile(path) && push!(files, path)
        end
    end
    entries = [
        string(relpath(path, root), ":", bytes2hex(sha256(read(path))))
            for path in sort!(files)
    ]
    push!(entries, string("Project.toml:", project_input_hash(project)))
    sort!(entries)
    _PACKAGE_INPUT_HASH[] = bytes2hex(sha256(join(entries, '\0')))
    return _PACKAGE_INPUT_HASH[]::String
end

"""Verify any locally available, ignored datasets against the committed
provenance declaration (`docs/notebooks/DATA_PROVENANCE.sha256`, one
`<sha256>  <path under docs/notebooks/data/>` per line). The file is optional,
and missing datasets are allowed: CI consumes the committed exports and does
not need private/large rendering inputs."""
function verify_local_data!()
    isfile(DATA_PROVENANCE) || return nothing
    for raw_line in eachline(DATA_PROVENANCE)
        line = strip(raw_line)
        (isempty(line) || startswith(line, '#')) && continue
        parts = split(line; limit = 2)
        length(parts) == 2 || error("malformed data provenance line: $(raw_line)")
        expected, relative = parts
        relative = strip(relative)
        path = joinpath(DATA_DIR, relative)
        isfile(path) || continue
        actual = bytes2hex(sha256(read(path)))
        actual == expected || error("dataset checksum mismatch: $(relative)")
    end
    return nothing
end

"Fingerprint both notebook input and the immutable exporter implementation."
function export_fingerprint(slug::AbstractString)
    snapshot_tree = string(Pkg.dependencies()[SNAPSHOT_UUID].tree_hash)
    build_lock = bytes2hex(sha256(read(joinpath(@__DIR__, "build_env", "Manifest.toml"))))
    docs_lock = bytes2hex(sha256(read(joinpath(@__DIR__, "Manifest.toml"))))
    driver = bytes2hex(sha256(read(@__FILE__)))
    payload = join(
        (
            source_hash(slug), package_input_hash(), snapshot_tree,
            build_lock, docs_lock, driver, EXPORT_CONTRACT,
        ), '\0'
    )
    return bytes2hex(sha256(payload))
end

"""Fingerprint recorded inside the cached HTML, or `nothing` if the cache is
missing or unreadable."""
function cached_hash(slug::AbstractString)::Union{String, Nothing}
    out = joinpath(STATIC_OUT, "$(slug).html")
    isfile(out) || return nothing
    html = read(out, String)
    m = match(
        Regex("<meta\\s+name=\"$(HASH_META)\"\\s+content=\"([0-9a-f]+)\""),
        html
    )
    return m === nothing ? nothing : m.captures[1]
end

needs_rebuild(slug::AbstractString)::Bool = cached_hash(slug) != export_fingerprint(slug)

"""Render one notebook → LEAN standalone HTML (Snapshot.jl therapy emit).

`export_notebook` opens the notebook in `session`, runs every cell (the
notebooks self-activate docs/Project.toml, so Metal/CairoMakie resolve as
always), SSRs the rendered outputs, attempts a WASM island per `@bind`
group (a group that can't compile falls back to its static output with a
visible note — never a blank), writes `<slug>.html` + `<slug>.islands/`,
and shuts the notebook down.

Returns the SHA-256 hash of the source after the run so the caller can diff
against the pre-render hash and report when Pluto's normalization shifted
the source.
"""
function render_notebook!(
        session::Pluto.ServerSession,
        src_path::AbstractString,
        dst_path::AbstractString
    )
    slug = splitext(basename(src_path))[1]
    html_path = export_notebook(
        src_path;
        output_dir = dirname(dst_path),
        therapy = true, # lean page: SSR cells + islands, NO statefile
        fragment = true, # native inline component for Therapy docs (no iframe)
        assets_base = "__VASCULATURE_NOTEBOOKS_BASE__",
        islands = true, # compile @bind groups to wasm where possible
        verify = true, # oracle-check compiled islands (node)
        optimize = :size,
        force_fallback_bonds = get(FORCE_FALLBACK_BONDS, slug, ()),
        env_dir = @__DIR__, # notebook imports resolve from docs/Project.toml
        session = session
    )
    post_hash = source_hash(slug)
    fingerprint = export_fingerprint(slug)
    # Embed the source+exporter fingerprint for the rebuild cache.
    meta = """<meta name="$(HASH_META)" content="$(fingerprint)">"""
    html = read(html_path, String)
    html = occursin("<head>", html) ? replace(html, "<head>" => "<head>" * meta; count = 1) :
        html * "\n" * meta * "\n"
    write(html_path, html)
    # export_notebook names the file after the notebook stem == our slug, so
    # html_path and dst_path coincide; assert so a future rename can't desync.
    abspath(html_path) == abspath(dst_path) ||
        error("export wrote $(html_path) but the gallery expects $(dst_path)")
    fragment_path = replace(html_path, r"\.html$" => ".fragment.html")
    isfile(fragment_path) || error("Snapshot did not emit expected fragment: $(fragment_path)")
    fragment = read(fragment_path, String)
    marker = "<!-- $(HASH_META): $(fingerprint) -->"
    write(fragment_path, marker * "\n" * fragment)
    return post_hash
end

"""Set of slugs the user wants to skip entirely (no render, not in the route
list). Comma-separated env var. Example:
    VASCULATURE_SKIP_NOTEBOOKS=00_cardiac_pipeline julia --project=docs docs/app.jl build
"""
function skipped_slugs()
    return Set(strip.(split(get(ENV, "VASCULATURE_SKIP_NOTEBOOKS", ""), ","; keepempty = false)))
end

"""Walk every notebook, rendering only those whose export fingerprint doesn't
match the cached fingerprint. Returns the full list of slugs (cached + freshly
rendered) so callers can register routes for all of them."""
function export_notebooks()
    verify_local_data!()
    slugs = notebook_slugs()
    skip = skipped_slugs()
    if !isempty(skip)
        kept = filter(s -> !(s in skip), slugs)
        dropped = filter(s -> s in skip, slugs)
        isempty(dropped) || println("[notebooks] VASCULATURE_SKIP_NOTEBOOKS → skipping $(join(dropped, ", "))")
        slugs = kept
    end
    isempty(slugs) && (println("[notebooks] no .jl files in $(NOTEBOOKS_DIR)"); return slugs)

    isdir(STATIC_OUT) || mkpath(STATIC_OUT)

    force = haskey(ENV, "VASCULATURE_FORCE_NB_REBUILD")
    to_render = force ? slugs : filter(needs_rebuild, slugs)

    if isempty(to_render)
        println("[notebooks] all $(length(slugs)) cached (export fingerprint matches) — no rebuild")
        return slugs
    end

    println(
        "[notebooks] rendering $(length(to_render))/$(length(slugs)) " *
            (force ? "(forced)" : "(source changed)") * "…"
    )

    # disable_writing_notebook_files=true: keep Pluto from re-saving the
    # source .jl during open(), so the hash taken before the render still
    # matches the file on disk afterwards.
    options = Pluto.Configuration.from_flat_kwargs(disable_writing_notebook_files = true)
    session = Pluto.ServerSession(; options = options)
    for slug in to_render
        src = joinpath(NOTEBOOKS_DIR, "$(slug).jl")
        dst = joinpath(STATIC_OUT, "$(slug).html")
        pre_hash = source_hash(slug)
        println("  ▸ $(slug)  (hash $(first(pre_hash, 12))…)")
        t0 = time()
        post_hash = render_notebook!(session, src, dst)
        if post_hash != pre_hash
            println("    ⚠ Pluto normalized source: hash $(first(pre_hash, 12))… → $(first(post_hash, 12))… (embedding post-render hash so cache stays consistent)")
        end
        kb = round(filesize(dst) / 1024; digits = 1)
        println("    ✓ $(round(time() - t0; digits = 1))s  →  $(kb) KB")
    end

    return slugs
end

if abspath(PROGRAM_FILE) == @__FILE__
    export_notebooks()
end
