# SPDX-License-Identifier: GPL-3.0-or-later
# Post-processing (src/postprocess.jl), a port of upstream's coronary_postprocess.py.
#
# Reference outputs, when present, come from running the upstream script on trees grown by
# this package (the directory is regenerable and not committed):
#   julia -t 8 --project -e 'using Vasculature; S = "oracle/build/src/samples";
#       generate("$S/case0.vol", "$S/case0.xml"; nterm = 2000, aperf = 10000, min_dist = 1,
#                root = (-26, 41, 13), seed = 1, T = Float64, export_xml = "oracle/out/pp/tree2000.xml")'
#   python -I oracle/pp_reference.py oracle/build/src/postprocessing/coronary_postprocess.py \
#       oracle/out/pp/tree2000.xml oracle/out/pp/ref2000 oracle/build/src/samples/case0.vol
# Every oracle/out/pp/tree<N>.xml with a matching ref<N>_*.csv is compared: the CSV files
# must be identical, and with ref<N>_internals.txt the per-segment orders and hemodynamics
# must agree bit for bit. The tube mesh is compared for the smallest reference tree. These
# comparisons switch on the upstream behaviour (`mass = :upstream, elements = :upstream`);
# the defaults correct two upstream bugs (see docs/postprocess.md).
using Vasculature
using Test

const PP_DIR = joinpath(@__DIR__, "..", "oracle", "out", "pp")
const PP_CSV = ("connectivity_matrix", "diameter_stats", "diameter_ranges", "hemodynamics")

function compare_internals(path, r)
    lines = readlines(path)
    rows = split.(lines[2:end])
    h = r.hemodynamics
    @test parse(Float64, split(lines[1])[2]) === h.mass_g
    @test [parse(Int, x[1]) for x in rows] == r.strahler
    @test [parse(Int, x[2]) + 1 for x in rows] == r.seg2elem
    @test [parse(Int, x[4]) for x in rows] == r.orders
    sol = (h.resting.flow, h.resting.pressure, h.resting.wss, h.stress.flow, h.stress.pressure, h.stress.wss)
    for k in 1:6
        @test [parse(Float64, x[4 + k]) for x in rows] == sol[k]
    end
    return
end

# vtkOBJWriter output vs ours: same structure; numbers equal up to the last bit of a Float32
# (VTK's arm64 build contracts a*b+c into FMAs, visible only in near-zero components).
function compare_obj(ref, ours)
    nl = ndiff = 0
    maxdiff = 0.0
    for (a, b) in zip(Iterators.drop(eachline(ref), 1), Iterators.drop(eachline(ours), 1))
        nl += 1
        a == b && continue
        ndiff += 1
        fa, fb = split(a), split(b)
        @test fa[1] == fb[1] && fa[1] in ("v", "vn")
        maxdiff = max(maxdiff, maximum(abs(parse(Float64, x) - parse(Float64, y)) for (x, y) in zip(fa[2:end], fb[2:end])))
    end
    @test nl == countlines(ref) - 1 == countlines(ours) - 1
    @test ndiff <= nl ÷ 100
    return @test maxdiff < 1.0e-6
end

@testset "post-processing" begin
    @testset "numpy and printf primitives" begin
        C = Vasculature
        @test C._np_median([3.0, 1.0, 2.0]) == 2.0
        @test C._np_median([4.0, 1.0, 2.0, 3.0]) == 2.5
        @test C._np_percentile([1.0, 2.0, 3.0, 4.0], 25) == 1.75
        @test C._np_percentile([1.0, 2.0, 3.0, 4.0], 75) == 3.25
        @test C._np_percentile([5.0], 25) == 5.0
        @test isnan(C._np_std([1.0]; ddof = 1))
        @test C._np_sum(fill(0.1, 1000)) == 100.00000000000001   # np.sum(np.full(1000, 0.1))
        @test C._pyfmt(NaN, 3) == "nan" && C._pyfmt(2.5, 0) == "2" && C._pyfmt(0.125, 2) == "0.12"
        @test [C._to_chars(x) for x in (0.0, 23.0, -26.52526092529297, 1.0e-5, 0.000123, 2.9943092158646323e-5, 1.0e21)] ==
            ["0", "23", "-26.52526092529297", "1e-05", "0.000123", "2.9943092158646323e-05", "1e+21"]
        @test C.literature_order_of(3.176) == 11 && C.literature_order_of(0.009) == 1
    end

    @testset "smoke test on a 200-terminal tree" begin
        E = grow_case0(200; T = Float64)
        dir = mktempdir()
        xml = export_gxl(E.m, joinpath(dir, "tree.xml"))
        prefix = joinpath(dir, "out", "pp")
        r = postprocess(xml, prefix; domain_vol = sample("case0.vol"), plots = false, verbose = false)
        n = nsegments(E.m) - 1                          # one XML edge per non-root-node segment
        @test length(r.orders) == length(r.diameters) == n
        @test all(1 .<= r.orders .<= 11)
        @test r.diameters == 2 .* parse.(Float64, Vasculature._g.(E.m.radius[2:end]))   # XML keeps %g
        for suffix in ("_tree.vtp", "_tree.obj", (string("_", k, ".csv") for k in PP_CSV)...)
            @test isfile(prefix * suffix)
        end
        @test readline(prefix * "_diameter_stats.csv") == "Order,Diameter (units) ± SD"
        vtp = read(prefix * "_tree.vtp", String)
        @test occursin("NumberOfLines=\"$n\"", vtp) && occursin("Name=\"final_order\"", vtp)
        @test countlines(prefix * "_tree.obj") == 1 + n * (4 * 24 * 2 + 2 * 24 + 2 * 22)
        # the network solve conserves flow, keeps 100 mmHg at the inlet side and drops pressure
        h = r.hemodynamics
        @test h.mass_g ≈ 254.98 * 1.05                  # 254,980 voxels > 128, 1 mm^3 each
        @test Vasculature.domain_mass_g(sample("case0.vol"); mass = :upstream)[1] ≈ 41.7396
        st = h.stress
        net = Vasculature._hemo_network(r.tree)
        @test st.flow[net.root] ≈ 3.5 * h.mass_g / 60
        for i in 1:n
            ch = net.children[i]
            isempty(ch) && continue
            @test sum(st.flow[ch]) ≈ st.flow[i]
            @test all(st.pressure[ch] .< st.pressure[i])
        end
        @test all(0 .< st.pressure .< 100)
        # default elements never hold two daughters of one bifurcation; upstream's can
        tp = r.tree
        node2seg = Dict(tp.to[i] => i for i in 1:n)
        par = [get(node2seg, tp.from[i], 0) for i in 1:n]
        siblings(o) = count(segs -> !allunique(par[segs]), o.elem_segs)
        @test siblings(Vasculature.diameter_defined_orders(tp)) == 0
        @test siblings(Vasculature.diameter_defined_orders(tp; elements = :upstream)) > 0
    end

    refs = isdir(PP_DIR) ? [m.captures[1] for m in filter(!isnothing, match.(r"^tree(\w+)\.xml$", readdir(PP_DIR)))] : String[]
    filter!(n -> all(isfile(joinpath(PP_DIR, "ref$(n)_$k.csv")) for k in PP_CSV), refs)
    if isempty(refs)
        @info "post-processing: no upstream reference outputs in oracle/out/pp; skipping the comparison"
    end
    smallest = isempty(refs) ? "" : refs[argmin(filesize.(joinpath.(PP_DIR, "tree" .* refs .* ".xml")))]
    @testset "matches upstream coronary_postprocess.py on tree$n" for n in refs
        dir = mktempdir()
        withobj = n == smallest && isfile(joinpath(PP_DIR, "ref$(n)_tree.obj"))
        r = postprocess(
            joinpath(PP_DIR, "tree$n.xml"), joinpath(dir, "jl"); obj = withobj,
            domain_vol = joinpath(PP_DIR, "..", "..", "build", "src", "samples", "case0.vol"),
            plots = false, verbose = false, mass = :upstream, elements = :upstream
        )
        for k in PP_CSV
            @test read(joinpath(dir, "jl_$k.csv"), String) == read(joinpath(PP_DIR, "ref$(n)_$k.csv"), String)
        end
        internals = joinpath(PP_DIR, "ref$(n)_internals.txt")
        isfile(internals) && compare_internals(internals, r)
        withobj && compare_obj(joinpath(PP_DIR, "ref$(n)_tree.obj"), joinpath(dir, "jl_tree.obj"))
    end
end
