# api_sections.jl — the section map of the generated API reference (docs/src/routes/api/index.jl).
#
# Pure data plus two functions, with no package dependencies. The docs route and the package's test suite
# (test/docs.jl) both include this file, so the page and the test agree on where every documented name goes.
#
# A name resolves to a section in two steps, first match wins:
#   1. `API_NAME_SECTION`: per-name overrides (the organ presets live in one file but belong with their stage);
#   2. `API_FILE_SECTION`: the path, relative to `src/`, of the file holding the name's docstring: an exact file, or
#      else the longest matching folder prefix.
# A name that matches neither lands in "Other" and fails test/docs.jl, so a new source file is mapped deliberately.

"""The documented modules, each with the prefix its names take on the page."""
const API_MODULES = [("Vasculature", ""), ("Flow", "Flow."), ("Body", "Body."), ("DynamicPhantom", "DynamicPhantom."), ("XCAT", "XCAT.")]

"""Ordered sections of the API page: `(id, title, intro)`. The intro is Markdown."""
const API_SECTIONS = [
    (
        id = "overview", title = "Overview", intro = """
        Every organ walks the same four stages, each one function and one spec struct:

        ```text
        GrowthSpec    ─▶ grow       ─▶ the grown trees
        FlowSpec      ─▶ flow       ─▶ every vessel's flow, every voxel's ground truth
        TransportSpec ─▶ transport  ─▶ the iodine in every vessel and bed, every second
        PhantomSpec   ─▶ phantom    ─▶ the CT phantom: labels, compartments, curves, ground truth
        ```

        An organ is one instance of each spec (`HEART_GROWTH`, `HEART_FLOW`, …). Everything below is generated from the
        package's docstrings when the site is built, so it cannot drift from the code, and each entry links to its source.
        The [Getting Started](../getting-started/) guide shows how the pieces fit together.
        """,
    ),
    (id = "grow", title = "Grow", intro = "Grow a vascular tree from the imaged inlets into the organ's tissue."),
    (id = "flow", title = "Flow", intro = "Blood flow in the grown bed from its two pressures, the vessel sizes and the arterioles' tone."),
    (id = "transport", title = "Transport", intro = "A contrast injection through the whole body and the organ's grown tree."),
    (id = "phantom", title = "Phantom", intro = "The voxelized dynamic phantom, and loading it for a CT simulator."),
    (id = "xcat", title = "XCAT Inputs", intro = "`Vasculature.XCAT` turns XCAT's surfaces, meshes and label volumes into the domains, inlet trees and blockers `grow` takes."),
    (id = "files", title = "Files and Checks", intro = "Tree and voxel files, morphometry and quality checks."),
    (id = "other", title = "Other", intro = "Names the section map does not place yet."),
]

"""The order names lead their sections in (the spec and its presets first); the rest follow by file and line."""
const API_LEADS = [
    "GrowthSpec", "HEART_GROWTH", "HEART_VEIN_GROWTH", "LUNG_GROWTH", "LUNG_VEIN_GROWTH", "grow",
    "FlowSpec", "HEART_FLOW", "LUNG_FLOW", "flow",
    "TransportSpec", "HEART_TRANSPORT", "LUNG_TRANSPORT", "transport",
    "PhantomSpec", "HEART_PHANTOM", "LUNG_PHANTOM", "phantom",
]

"""Per-name overrides, checked first."""
const API_NAME_SECTION = Dict{String, String}(
    "HEART_GROWTH" => "grow", "HEART_VEIN_GROWTH" => "grow", "LUNG_GROWTH" => "grow", "LUNG_VEIN_GROWTH" => "grow",
    "HEART_TRANSPORT" => "transport", "LUNG_TRANSPORT" => "transport",
    "HEART_PHANTOM" => "phantom", "LUNG_PHANTOM" => "phantom",
    "XCAT_BODY_POOLS" => "phantom", "XCAT_POOL_OF_LABEL" => "phantom", "XCAT_STATIC" => "phantom",
)

"""Source file or folder (relative to `src/`) → section, checked second (an exact file, or the longest folder)."""
const API_FILE_SECTION = Dict{String, String}(
    "spec.jl" => "grow", "grow.jl" => "grow", "model.jl" => "grow", "params.jl" => "grow", "domain.jl" => "grow",
    "view.jl" => "grow", "generate.jl" => "grow", "batched.jl" => "grow",
    "perfusion/flow.jl" => "flow", "perfusion/pipeline.jl" => "flow", "perfusion/organs.jl" => "flow",
    "perfusion/territory.jl" => "flow",
    "perfusion/transport.jl" => "transport", "perfusion/body.jl" => "transport",
    "phantom.jl" => "phantom", "perfusion/dynamic_phantom.jl" => "phantom", "perfusion/organct.jl" => "phantom",
    "xcat/" => "xcat",
    "io_gxl.jl" => "files", "io_vol.jl" => "files", "qc.jl" => "files", "morphometry.jl" => "files",
    "postprocess.jl" => "files", "postprocess/" => "files",
)

"""
    api_section(name, src_path) -> String

Section id for documented `name` (without its module prefix), whose docstring lives in `src_path` (relative to `src/`,
`/` separators; `nothing` when unknown). Returns `"other"` when nothing matches.
"""
function api_section(name::AbstractString, src_path::Union{Nothing, AbstractString})
    haskey(API_NAME_SECTION, name) && return API_NAME_SECTION[name]
    src_path === nothing && return "other"
    haskey(API_FILE_SECTION, src_path) && return API_FILE_SECTION[src_path]
    best = ""
    for p in keys(API_FILE_SECTION)
        endswith(p, "/") && startswith(src_path, p) && length(p) > length(best) && (best = p)
    end
    return isempty(best) ? "other" : API_FILE_SECTION[best]
end

"""
    api_doc_location(docstr_data) -> (path, line)

`(path relative to src/, line)` of a `Base.Docs.DocStr`'s `data` dict, or `(nothing, 0)`.
"""
function api_doc_location(data)
    p = get(data, :path, nothing)
    p === nothing && return (nothing, 0)
    s = replace(string(p), '\\' => '/')
    i = findlast("/src/", s)
    rel = i === nothing ? basename(s) : s[(last(i) + 1):end]
    return (rel, Int(get(data, :linenumber, 0)))
end
