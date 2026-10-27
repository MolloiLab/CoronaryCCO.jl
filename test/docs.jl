# Docs coverage: every public name has a docstring, and the API reference's section map (docs/src/api_sections.jl, the
# file the generated /api/ page uses) places every one of them. Dependency-free: the map is plain Julia data plus two
# functions.

const _APIMap = Module(:_APIMap)
Base.include(_APIMap, joinpath(pkgdir(Vasculature), "docs", "src", "api_sections.jl"))

# the names the API page documents, as it collects them: each module's exports (XCAT: its documented public names),
# a submodule's name documented once where the package also exports it
function _documented_names()
    top = Vasculature
    out = Tuple{Module, Symbol}[]
    for (m, prefix) in _APIMap.API_MODULES
        M = m == "Vasculature" ? top : getfield(top, Symbol(m))
        ns = prefix == "XCAT." ?
            [b.var for b in keys(Base.Docs.meta(M)) if b.mod === M && !startswith(string(b.var), "_") && isdefined(M, b.var)] :
            filter(n -> n !== nameof(M), names(M))
        M === top || filter!(n -> !(n in names(top) && isdefined(top, n) && getfield(top, n) === getfield(M, n)), ns)
        append!(out, [(M, n) for n in unique(ns)])
    end
    return out
end

@testset "docs: every public name has a docstring and a section" begin
    undocumented = String[]
    unplaced = String[]
    for (M, n) in _documented_names()
        b = Base.Docs.Binding(M, n)
        meta = Base.Docs.meta(b.mod)
        if !haskey(meta, b)
            push!(undocumented, "$(nameof(M)).$n")
            continue
        end
        md = meta[b]
        path, _ = _APIMap.api_doc_location(md.docs[first(md.order)].data)
        _APIMap.api_section(string(n), path) == "other" && push!(unplaced, "$(nameof(M)).$n ($(something(path, "?")))")
    end
    @test undocumented == String[]
    @test unplaced == String[]
end
