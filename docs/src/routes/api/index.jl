# /api/ — the API reference, GENERATED at build time from Vasculature's docstrings.
#
# Nothing on this page is hand-written about the API: every exported name (and every documented name of the XCAT, Flow,
# Body and DynamicPhantom modules) is listed with its signature, its rendered docstring and a link to its source. The
# section a name appears in comes from docs/src/api_sections.jl, shared with test/docs.jl; a name the map does not place
# lands in "Other" and fails the test suite. The name filter is the ApiFilter island (src/components/ApiFilter.jl).
#
# Therapy includes this file into its `TherapyApp` module; its last expression is the route.

import Markdown
import Vasculature

include(joinpath(@__DIR__, "..", "..", "api_sections.jl"))

# The package is not registered or tagged, so source links point at `main`.
const _API_GITHUB_SRC = "https://github.com/MolloiLab/Vasculature.jl/blob/main/src"

_api_esc(s) = replace(string(s), '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;")
_api_unesc(s) = replace(
    replace(s, r"&#(\d+);" => m -> string(Char(parse(Int, m[3:(end - 1)])))),
    "&lt;" => "<", "&gt;" => ">", "&quot;" => "\"", "&amp;" => "&",
)

# `[`name`](@ref)` → an in-page link when `name` is on this page, plain code otherwise.
function _api_ref_target(txt, onpage)
    t = replace(_api_unesc(txt), r"^(Vasculature\.)" => "")
    t = strip(first(split(t, r"[\(\{\s]"; limit = 2)))
    return t in onpage ? String(t) : nothing
end

function _api_style_html(html::AbstractString, onpage::Set{String})
    h = replace(
        html, r"<a href=\"@ref[^\"]*\"><code>(.*?)</code></a>" => function (m)
            inner = match(r"<code>(.*?)</code>", m).captures[1]
            tgt = _api_ref_target(inner, onpage)
            return tgt === nothing ? "<code>$inner</code>" : "<a href=\"#$(_api_esc(tgt))\"><code>$inner</code></a>"
        end
    )
    h = replace(h, r"<a href=\"@ref[^\"]*\">(.*?)</a>" => s"\1")
    h = replace(h, "<pre><code>" => "<pre><code class=\"language-text\">")
    h = replace(h, r"<a href=\"(?!#)" => "<a target=\"_blank\" rel=\"noopener\" href=\"")
    return replace(h, r"<h[1-6]>" => "<h4>", r"</h[1-6]>" => "</h4>")
end

_api_md_html(md, onpage) = _api_style_html(Markdown.html(md), onpage)

# ── collecting the docstrings ─────────────────────────────────────────────────────────────────────────────────────────

function _api_kind(M::Module, n::Symbol)
    isdefined(M, n) || return "binding"
    v = getfield(M, n)
    v isa Module && return "module"
    if v isa Type
        t = Base.unwrap_unionall(v)
        t isa DataType || return "type"
        isabstracttype(t) && return "abstract type"
        return ismutabletype(t) ? "mutable struct" : "struct"
    end
    v isa Function && return startswith(string(n), "@") ? "macro" : "function"
    return "constant"
end

# Julia Markdown reads an underscore inside a word (`tissue_ml`, `r_root`) as emphasis, so prose written without
# backticks would lose its underscores. Escape every intra-word underscore outside code.
function _api_escape_underscores(text::AbstractString)
    out = IOBuffer()
    fenced = false
    for line in split(text, '\n'; keepempty = true)
        if startswith(lstrip(line), "```")
            fenced = !fenced
            println(out, line)
            continue
        end
        if fenced || startswith(line, "    ") || startswith(line, "\t")
            println(out, line)
            continue
        end
        parts = split(line, r"(`+)"; keepempty = true)
        ticks = [m.match for m in eachmatch(r"`+", line)]
        incode = false
        opener = ""
        for (i, part) in enumerate(parts)
            print(out, incode ? part : replace(part, r"(?<=[\p{L}\p{N}])_(?=[\p{L}\p{N}])" => "\\_"))
            i <= length(ticks) || break
            t = ticks[i]
            if !incode
                incode = true
                opener = t
            elseif t == opener
                incode = false
            end
            print(out, t)
        end
        println(out)
    end
    return String(chomp(String(take!(out))))
end

_api_parse(d) = isempty(d.text) ? d.object : Markdown.parse(_api_escape_underscores(join(string.(d.text))))
_api_blocks(md) = length(md.content) == 1 && md.content[1] isa Markdown.MD ? _api_blocks(md.content[1]) : md.content

function _api_fallback_signature(M, n, kind)
    kind == "module" && return "module $(nameof(getfield(M, n)))"
    kind in ("struct", "mutable struct", "abstract type", "type") && return "$(kind) $(n)"
    kind == "constant" && return "$(n)::$(typeof(getfield(M, n)))"
    if kind == "function"
        ms = [m for m in methods(getfield(M, n)) if parentmodule(m) === M]
        sigs = unique([replace(sprint(show, m), r"\s+@\s.*$"s => "", r"Vasculature\.(\w+\.)?" => "") for m in ms])
        isempty(sigs) || return join(sigs, "\n")
    end
    return string(n)
end

# The names to document from one module: its exports, or (XCAT, which exports nothing) every documented public name.
function _api_names(M::Module, prefix::AbstractString, top::Module)
    meta = Base.Docs.meta(M)
    if prefix == "XCAT."
        ns = [b.var for b in keys(meta) if b.mod === M && !startswith(string(b.var), "_") && isdefined(M, b.var)]
    else
        ns = filter(n -> n !== nameof(M), names(M))
    end
    # a submodule's name that the package also exports at top level is documented once, there
    M === top || filter!(n -> !(isdefined(top, n) && n in names(top) && getfield(top, n) === getfield(M, n)), ns)
    return sort!(unique(ns); by = n -> lowercase(string(n)))
end

"""Every documented name: its section, kind, and one `(sig, body, path, line)` per docstring."""
function _api_collect()
    top = Vasculature
    modules = [(m == "Vasculature" ? top : getfield(top, Symbol(m)), p) for (m, p) in API_MODULES]
    listed = [(M, p, n) for (M, p) in modules for n in _api_names(M, p, top)]
    onpage = Set(string(p, n) for (_, p, n) in listed)
    entries = []
    for (M, prefix, n) in listed
        kind = _api_kind(M, n)
        b = Base.Docs.Binding(M, n)
        meta = Base.Docs.meta(b.mod)
        docs = []
        if haskey(meta, b)
            md = meta[b]
            for sig in md.order
                d = md.docs[sig]
                blocks = copy(_api_blocks(_api_parse(d)))
                signature = nothing
                if !isempty(blocks) && blocks[1] isa Markdown.Code && blocks[1].language in ("", "julia") &&
                        occursin(replace(string(n), "!" => ""), blocks[1].code)
                    signature = blocks[1].code
                    popfirst!(blocks)
                end
                path, line = api_doc_location(d.data)
                push!(docs, (sig = signature, body = _api_md_html(Markdown.MD(blocks), onpage), path = path, line = line))
            end
        end
        path = isempty(docs) ? nothing : docs[1].path
        name = string(prefix, n)
        push!(
            entries, (
                name = name, kind = kind, section = api_section(string(n), path), docs = docs, path = path,
                line = isempty(docs) ? 0 : docs[1].line, fallback = _api_fallback_signature(M, n, kind),
            )
        )
    end
    return entries, onpage
end

# ── rendering ────────────────────────────────────────────────────────────────────────────────────────────────────────

_api_srcurl(path, line) = "$(_API_GITHUB_SRC)/$(path)" * (line > 0 ? "#L$(line)" : "")

function _api_source_link(path, line)
    path === nothing && return ""
    return """<a href="$(_api_esc(_api_srcurl(path, line)))" target="_blank" rel="noopener" title="src/$(_api_esc(path))$(line > 0 ? ":$line" : "") on GitHub" class="tv-api-src">source ↗</a>"""
end

_api_sig_html(code) = """<pre class="tv-api-sig"><code class="language-julia">$(_api_esc(code))</code></pre>"""

function _api_entry_html(e)
    io = IOBuffer()
    id = _api_esc(e.name)
    print(io, """<article id="$id" data-api-entry="$(_api_esc(lowercase(e.name)))" class="tv-api-entry">""")
    print(io, """<header class="tv-api-head"><a href="#$id" class="tv-api-name">$id</a><span class="tv-api-kind">$(_api_esc(e.kind))</span><span class="flex-1"></span>""")
    length(e.docs) == 1 && print(io, _api_source_link(e.path, e.line))
    print(io, "</header>")
    isempty(e.docs) && print(io, _api_sig_html(e.fallback), """<p class="tv-api-missing">No docstring yet.</p>""")
    for (i, d) in enumerate(e.docs)
        i > 1 && print(io, """<hr class="tv-api-rule">""")
        length(e.docs) > 1 && print(io, """<div class="flex justify-end">""", _api_source_link(d.path, d.line), "</div>")
        print(io, _api_sig_html(something(d.sig, e.fallback)))
        print(io, """<div class="tv-api-doc">""", d.body, "</div>")
    end
    print(io, "</article>")
    return String(take!(io))
end

function _api_section_html(sec, entries, onpage)
    io = IOBuffer()
    print(io, """<div class="space-y-4"><div class="tv-api-doc tv-api-intro">""", _api_md_html(Markdown.parse(sec.intro), onpage), "</div>")
    if !isempty(entries)
        print(io, """<div class="flex flex-wrap gap-1.5 pt-1">""")
        for e in entries
            print(io, """<a href="#$(_api_esc(e.name))" data-api-chip="$(_api_esc(lowercase(e.name)))" class="tv-api-chip">$(_api_esc(e.name))</a>""")
        end
        print(io, "</div>")
        lastpath = missing
        for e in entries
            if !isequal(e.path, lastpath)
                lastpath = e.path
                print(io, """<p data-api-file class="tv-api-file">$(_api_esc(e.path === nothing ? "undocumented" : "src/" * e.path))</p>""")
            end
            print(io, _api_entry_html(e))
        end
    end
    print(io, "</div>")
    return String(take!(io))
end

# The name filter (the pattern of BasisSimulator.jl's docs): a small dependency-free script. `data-therapy-rerun` makes
# Therapy's client router run it again after a navigation swaps this page in.
const _API_FILTER_JS = """
<script data-therapy-rerun>(function(){var i=document.getElementById('api-filter');if(!i||i.dataset.bound)return;i.dataset.bound='1';
var es=[].slice.call(document.querySelectorAll('[data-api-entry]')),cs=[].slice.call(document.querySelectorAll('[data-api-chip]')),
ss=[].slice.call(document.querySelectorAll('[data-api-section]')),fs=[].slice.call(document.querySelectorAll('[data-api-file]')),c=document.getElementById('api-filter-count');
function run(){var q=i.value.trim().toLowerCase(),n=0;
es.forEach(function(e){var m=!q||e.getAttribute('data-api-entry').indexOf(q)>=0;e.style.display=m?'':'none';if(m)n++;});
cs.forEach(function(e){e.style.display=(!q||e.getAttribute('data-api-chip').indexOf(q)>=0)?'':'none';});
fs.forEach(function(f){f.style.display=q?'none':'';});
ss.forEach(function(s){var any=false;s.querySelectorAll('[data-api-entry]').forEach(function(e){if(e.style.display!=='none')any=true;});s.style.display=(q&&!any)?'none':'';});
if(c)c.textContent=q?(n+' of '+es.length+' names match'):(es.length+' names');}
i.addEventListener('input',run);
document.addEventListener('keydown',function(ev){var a=document.activeElement;if(ev.key==='/'&&a!==i&&!/^(input|textarea|select)\$/i.test((a&&a.tagName)||'')){ev.preventDefault();i.focus();}});
run();})();</script>
"""

let BASE = get(ENV, "VASCULATURE_BASE", "")
    () -> begin
        entries, onpage = _api_collect()
        by_section = Dict{String, Vector{Any}}()
        for e in entries
            push!(get!(by_section, e.section, Any[]), e)
        end
        lead_rank(n) = (i = findfirst(==(n), API_LEADS); i === nothing ? typemax(Int) : i)
        for v in values(by_section)
            sort!(v; by = e -> (lead_rank(e.name), something(e.path, "~"), e.line, e.name))
        end
        shown = [s for s in API_SECTIONS if s.id == "overview" || !isempty(get(by_section, s.id, []))]
        toc = [(s.id, s.title) for s in shown]

        PageWithTOC(
            toc, Div(
                :class => "max-w-3xl space-y-10",
                Div(
                    :class => "space-y-4",
                    H1(:class => "tv-display text-5xl sm:text-6xl", "API Reference"),
                    P(
                        :class => "tv-body text-[17px]",
                        "Every public name of Vasculature.jl, generated from its docstrings when the site is built. ",
                        "Four stages, each one function and one spec struct, and the tools around them. The ",
                        A(:href => "$(BASE)/getting-started/", :class => "tv-link", "Getting Started"),
                        " guide and the ",
                        A(:href => "$(BASE)/examples/", :class => "tv-link", "Examples"),
                        " show them in real pipelines."
                    ),
                ),
                # the filter bar is a direct child of the page column, so it stays pinned under the header
                Div(
                    :class => "tv-search-bar",
                    Div(
                        :class => "flex items-center gap-3",
                        Input(
                            :id => "api-filter", :type => "search", :class => "tv-search",
                            :placeholder => "Filter names…  (press / to focus)", :autocomplete => "off", :spellcheck => "false",
                            Symbol("aria-label") => "Filter the API reference by name",
                        ),
                        Span(:id => "api-filter-count", :class => "shrink-0 text-xs font-mono tv-muted", "$(length(entries)) names"),
                    ),
                ),
                [
                    Section(
                        :id => s.id, :class => "space-y-4 scroll-mt-24", Symbol("data-api-section") => s.id,
                        H2(:class => "tv-title text-3xl", s.title),
                        RawHtml(_api_section_html(s, get(by_section, s.id, []), onpage)),
                    ) for s in shown
                ]...,
                RawHtml(_API_FILTER_JS),
            )
        )
    end
end
