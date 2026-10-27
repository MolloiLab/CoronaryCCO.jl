#!/usr/bin/env julia
# Reject broken, stale, or machine-specific committed notebook exports.
#
# Usage (stdlib only, no environment needed):
#   julia docs/verify_notebook_exports.jl

using SHA: bytes2hex, sha256
using TOML

const DOCS = @__DIR__
const ROOT = joinpath(DOCS, "notebooks-static")
const SOURCES = joinpath(DOCS, "notebooks")
const DATA_PROVENANCE = joinpath(SOURCES, "DATA_PROVENANCE.sha256")
const HASH_META = "vasculature-export-fingerprint"
# Must match EXPORT_CONTRACT in extract_all.jl byte for byte.
const EXPORT_CONTRACT = "vasculature-lean-v1|therapy=true|fragment=true|islands=true|verify=true|optimize=size|forced-fallback="
const BAD_HTML = ["rendered Pluto error" => "<jlerror", "serialized Pluto exception" => "plain_error"]
const BAD_TEXT = ["developer-machine path" => "/users/daleblack/", "CI workspace path" => "/home/runner/work/"]
const BAD_TEXT_PATTERNS = [
    "macOS user path" => r"/(?:private/)?users/[^/\s<]+/",
    "Unix home path" => r"/home/[^/\s<]+/",
    "temporary build path" => r"/(?:private/)?tmp/[^\s<\"']+",
    "mounted volume path" => r"/volumes/[^\s<\"']+",
]
const BAD_REPORT_REASONS = ("failed to parse", "package ", "not found in current path")
const TEXT_SUFFIXES = (".html", ".json", ".js", ".css", ".map", ".txt")
# slug => (bond, Set(cell ids)): mirrors FORCE_FALLBACK_BONDS in extract_all.jl.
const FORCED_FALLBACKS = Dict{String, Tuple{String, Set{String}}}()

source_hash(path) = bytes2hex(sha256(read(path)))

"Hash Project.toml while ignoring the package's release-only version field."
function project_input_hash(path)
    normalized = replace(read(path, String), r"(?m)^version\s*=\s*\"[^\"]*\"\s*$" => "version = \"<release-version>\"")
    return bytes2hex(sha256(normalized))
end

"Digest committed package code and the (optional) committed data provenance."
function package_input_hash()
    root = dirname(DOCS)
    files = isfile(DATA_PROVENANCE) ? [DATA_PROVENANCE] : String[]
    for (base, _, names) in walkdir(joinpath(root, "src")), name in names
        path = joinpath(base, name)
        isfile(path) && push!(files, path)
    end
    entries = [string(relpath(path, root), ":", source_hash(path)) for path in files]
    push!(entries, string("Project.toml:", project_input_hash(joinpath(root, "Project.toml"))))
    return bytes2hex(sha256(join(sort!(entries), '\0')))
end

"Snapshot's locked source tree in docs/build_env/Manifest.toml."
function snapshot_tree_hash()
    manifest = TOML.parsefile(joinpath(DOCS, "build_env", "Manifest.toml"))
    return only(manifest["deps"]["Snapshot"])["git-tree-sha1"]
end

function export_fingerprint(path, snapshot_tree, package_hash)
    payload = join(
        (
            source_hash(path), package_hash, snapshot_tree,
            source_hash(joinpath(DOCS, "build_env", "Manifest.toml")),
            source_hash(joinpath(DOCS, "Manifest.toml")),
            source_hash(joinpath(DOCS, "extract_all.jl")), EXPORT_CONTRACT,
        ), '\0'
    )
    return bytes2hex(sha256(payload))
end

function recorded_hash(path, fragment::Bool)
    pattern = fragment ? Regex("<!--\\s*$(HASH_META):\\s*([0-9a-f]{64})\\s*-->") :
        Regex("<meta\\s+name=\"$(HASH_META)\"\\s+content=\"([0-9a-f]{64})\"")
    m = match(pattern, read(path, String))
    return m === nothing ? nothing : m.captures[1]
end

# ── minimal JSON reader (objects → Dict, arrays → Vector) for the island reports ──
function json_parse(s::AbstractString)
    v, i = _json_value(s, _json_ws(s, firstindex(s)))
    _json_ws(s, i) > lastindex(s) || error("trailing JSON content")
    return v
end
_json_ws(s, i) = (
    while i <= lastindex(s) && isspace(s[i])
        i = nextind(s, i)
    end; i
)
function _json_value(s, i)
    c = s[i]
    if c == '{'
        d = Dict{String, Any}(); i = _json_ws(s, nextind(s, i))
        s[i] == '}' && return d, nextind(s, i)
        while true
            k, i = _json_string(s, _json_ws(s, i)); i = _json_ws(s, i)
            s[i] == ':' || error("expected ':' in JSON")
            v, i = _json_value(s, _json_ws(s, nextind(s, i))); d[k] = v; i = _json_ws(s, i)
            s[i] == '}' && return d, nextind(s, i)
            s[i] == ',' || error("expected ',' in JSON object")
            i = nextind(s, i)
        end
    elseif c == '['
        a = Any[]; i = _json_ws(s, nextind(s, i))
        s[i] == ']' && return a, nextind(s, i)
        while true
            v, i = _json_value(s, _json_ws(s, i)); push!(a, v); i = _json_ws(s, i)
            s[i] == ']' && return a, nextind(s, i)
            s[i] == ',' || error("expected ',' in JSON array")
            i = nextind(s, i)
        end
    elseif c == '"'
        return _json_string(s, i)
    elseif startswith(SubString(s, i), "true")
        return true, i + 4
    elseif startswith(SubString(s, i), "false")
        return false, i + 5
    elseif startswith(SubString(s, i), "null")
        return nothing, i + 4
    else
        m = match(r"-?\d+(\.\d+)?([eE][+-]?\d+)?", s, i)
        (m === nothing || m.offset != i) && error("invalid JSON value")
        n = m.captures[1] === nothing && m.captures[2] === nothing ? parse(Int, m.match) : parse(Float64, m.match)
        return n, i + ncodeunits(m.match)
    end
end
function _json_string(s, i)
    s[i] == '"' || error("expected JSON string")
    io = IOBuffer(); i = nextind(s, i)
    while s[i] != '"'
        if s[i] == '\\'
            i = nextind(s, i); e = s[i]
            if e == 'u'
                write(io, Char(parse(UInt16, s[(i + 1):(i + 4)]; base = 16))); i += 4
            else
                write(io, get(Dict('n' => '\n', 't' => '\t', 'r' => '\r', 'b' => '\b', 'f' => '\f'), e, e))
            end
        else
            write(io, s[i])
        end
        i = nextind(s, i)
    end
    return String(take!(io)), nextind(s, i)
end
read_json(path) = json_parse(read(path, String))

"A configured fallback must publish exactly one disabled group."
function check_forced_fallback(assets, bond, expected_cells)
    failures = String[]
    rel = relpath(assets, ROOT)
    missing_files = sort([n for n in ("report.json", "coverage.json", "islands.json", "shim.js") if !isfile(joinpath(assets, n))])
    for n in missing_files
        push!(failures, "$(rel)/$(n): missing configured-fallback asset")
    end
    isempty(missing_files) || return failures
    reason = "configured fallback (matched @bind $(bond)); island inference intentionally skipped"
    groups, manifest, coverage = try
        read_json(joinpath(assets, "report.json")), read_json(joinpath(assets, "islands.json")),
            read_json(joinpath(assets, "coverage.json"))
    catch
        return ["$(rel): invalid configured-fallback JSON"]
    end
    if length(groups) != 1
        push!(failures, "$(rel)/report.json: expected one group")
    else
        group = groups[1]
        cells = get(group, "cells", Any[])
        if get(group, "bonds", nothing) != [bond] || get(group, "judgement", nothing) != "fallback" ||
                get(group, "fallback_kind", nothing) != "configured" || get(group, "reasons", nothing) != [reason] ||
                Set(string(get(c, "id", nothing)) for c in cells) != expected_cells ||
                any(get(c, "ok", nothing) !== false for c in cells) ||
                any(get(c, "reasons", nothing) != [reason] for c in cells)
            push!(failures, "$(rel)/report.json: configured fallback contract mismatch")
        end
    end
    expected_runtime = [Dict("bonds" => [bond], "judgement" => "fallback", "fallback_kind" => "configured")]
    if get(manifest, "groups", nothing) != [] || get(manifest, "fallback_groups", nothing) != expected_runtime
        push!(failures, "$(rel)/islands.json: runtime fallback index mismatch")
    end
    n = length(expected_cells)
    expected_coverage = Dict(
        "groups" => Dict("island" => 0, "partial" => 0, "fallback" => 1, "total" => 1),
        "cells" => Dict("interactive" => 0, "fallback" => n, "total" => n),
    )
    coverage == expected_coverage || push!(failures, "$(rel)/coverage.json: configured fallback coverage mismatch")
    return failures
end

root_html() = isdir(ROOT) ? [joinpath(ROOT, f) for f in readdir(ROOT) if endswith(f, ".html")] : String[]

function main()
    sources = isdir(SOURCES) ? sort([joinpath(SOURCES, f) for f in readdir(SOURCES) if endswith(f, ".jl")]) : String[]
    if isempty(sources)
        # No notebooks yet: the site builds with an empty Examples page.
        stray = root_html()
        if !isempty(stray)
            println(stderr, "orphan exports without source notebooks:")
            foreach(p -> println(stderr, "  - ", basename(p)), sort(stray))
            return 1
        end
        println("no notebook sources in $(SOURCES); nothing to verify")
        return 0
    end

    failures = String[]
    expected_pages = Set{String}()
    snapshot_tree = try
        snapshot_tree_hash()
    catch
        println(stderr, "invalid Snapshot export lock in docs/build_env/Manifest.toml")
        return 1
    end
    package_hash = package_input_hash()
    for source in sources
        slug = splitext(basename(source))[1]
        digest = export_fingerprint(source, snapshot_tree, package_hash)
        for (suffix, fragment) in ((".html", false), (".fragment.html", true))
            page = joinpath(ROOT, slug * suffix)
            push!(expected_pages, page)
            if !isfile(page)
                push!(failures, "$(basename(page)): missing export")
                continue
            end
            found = recorded_hash(page, fragment)
            found == digest || push!(
                failures,
                "$(basename(page)): stale export (source $(first(digest, 12)), export $(found === nothing ? "none" : first(found, 12)))"
            )
            text = lowercase(read(page, String))
            for (label, needle) in BAD_HTML
                occursin(lowercase(needle), text) && push!(failures, "$(basename(page)): $(label) ($(repr(needle)))")
            end
            fragment && occursin("localstorage.getitem('snap-theme')", text) &&
                push!(failures, "$(basename(page)): embedded fragment overrides host theme")
        end

        assets = joinpath(ROOT, "$(slug).islands")
        if haskey(FORCED_FALLBACKS, slug)
            bond, cells = FORCED_FALLBACKS[slug]
            append!(failures, check_forced_fallback(assets, bond, cells))
        end

        report = joinpath(assets, "report.json")
        if isfile(report)
            groups = try
                read_json(report)
            catch
                push!(failures, "$(relpath(report, ROOT)): invalid JSON")
                continue
            end
            for group in groups
                reasons = collect(Any, get(group, "reasons", Any[]))
                for cell in get(group, "cells", Any[])
                    append!(reasons, get(cell, "reasons", Any[]))
                end
                for reason in reasons
                    low = lowercase(string(reason))
                    any(n -> occursin(n, low), BAD_REPORT_REASONS) &&
                        push!(failures, "$(relpath(report, ROOT)): invalid island diagnostic")
                end
            end
        end
    end

    for orphan in sort(setdiff(root_html(), expected_pages))
        push!(failures, "$(basename(orphan)): orphan export without source notebook")
    end

    # Every textual artifact under notebooks-static is committed and copied by
    # Therapy. Catch machine paths even outside HTML (reports, maps, shims).
    if isdir(ROOT)
        for (base, _, names) in walkdir(ROOT), name in sort(names)
            any(s -> endswith(lowercase(name), s), TEXT_SUFFIXES) || continue
            artifact = joinpath(base, name)
            text = lowercase(read(artifact, String))
            for (label, needle) in BAD_TEXT
                occursin(needle, text) && push!(failures, "$(relpath(artifact, ROOT)): $(label) ($(repr(needle)))")
            end
            for (label, pattern) in BAD_TEXT_PATTERNS
                m = match(pattern, text)
                m === nothing || push!(failures, "$(relpath(artifact, ROOT)): $(label) ($(repr(m.match)))")
            end
        end
    end

    if !isempty(failures)
        println(stderr, "broken committed notebook exports:")
        foreach(f -> println(stderr, "  - ", f), failures)
        return 1
    end
    println("verified $(length(sources)) notebook source/export pairs")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
