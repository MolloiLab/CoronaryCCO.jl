let BASE = get(ENV, "VASCULATURE_BASE", "")
    () -> begin
        code_block = "tv-pre"
        prose = "tv-body text-[17px]"
        h2 = "tv-title text-3xl pt-6"
        inline = "tv-code"

        sections = [
            ("install", "Install"),
            ("pipeline", "The Four Stages in Code"),
            ("load", "Load the Phantom"),
            ("whats-next", "What's Next"),
        ]

        PageWithTOC(
            sections, Div(
                :class => "max-w-3xl space-y-6",
                H1(
                    :class => "tv-display text-5xl sm:text-6xl pb-2",
                    "Getting Started"
                ),
                P(
                    :class => prose,
                    "Vasculature.jl runs on Julia 1.12 or later. Growth runs on the CPU or on any GPU backend ",
                    "(Metal for Apple silicon, CUDA for NVIDIA, AMDGPU for AMD)."
                ),

                H2(:id => "install", :class => h2, "Install"),
                P(:class => prose, "From the Julia REPL:"),
                Pre(
                    :class => code_block,
                    Code(
                        :class => "language-julia", """using Pkg
                        Pkg.add(url = "https://github.com/MolloiLab/Vasculature.jl")"""
                    )
                ),
                P(:class => prose, "Then add your GPU backend, or skip it to grow on the CPU:"),
                Pre(
                    :class => code_block,
                    Code(
                        :class => "language-julia", """Pkg.add("Metal") # Apple silicon
                        Pkg.add("CUDA") # NVIDIA
                        Pkg.add("AMDGPU") # AMD"""
                    )
                ),

                H2(:id => "pipeline", :class => h2, "The Four Stages in Code"),
                P(
                    :class => prose,
                    "A sketch of the heart from organ to phantom. The organ surface (",
                    Code(:class => inline, "vertices"), ", ", Code(:class => inline, "triangles"),
                    "), the chambers to avoid, the tissue voxel centres and the voxel classifier come from your own ",
                    "segmentation. The lungs use ", Code(:class => inline, "LUNG_GROWTH"), ", ",
                    Code(:class => inline, "LUNG_FLOW"), ", ", Code(:class => inline, "LUNG_TRANSPORT"), " and ",
                    Code(:class => inline, "LUNG_PHANTOM"), " in the same calls."
                ),
                Pre(
                    :class => code_block,
                    Code(
                        :class => "language-julia", """using Vasculature
                        using Metal # or CUDA, AMDGPU; leave out for the CPU

                        # 1. Grow: organ surface → 1 mm domain, then a tree from each imaged inlet
                        domain = mesh_domain("heart.vol", vertices, triangles)
                        spec = GrowthSpec(HEART_GROWTH; device = MtlArray)
                        avoid = (chambers, origin_mm, 0.4) # forbidden mask, centre of its first cell (mm), cell size (mm)
                        tree = grow(spec, domain, ["lca.xml", "rca.xml"]; blockers = avoid)
                        veintree = grow(GrowthSpec(HEART_VEIN_GROWTH; device = MtlArray), domain, ["veins.xml"]) # veins grow from their imaged outlet trees

                        # 2. Flow: from the pressures, the vessel sizes and the arterioles' tone
                        arteries, veins = vtree(tree), vtree(veintree)
                        F = flow(HEART_FLOW, arteries, veins, tissue) # tissue: 3 × n voxel centres (mm)
                        F.rest.resist.mbf, F.stress.resist.mbf # mL/min/g at rest and at maximal dilation

                        # 3. Transport: a clinical injection through the body and the grown tree
                        Trest = transport(HEART_TRANSPORT, F.rest; times = HEART_PHANTOM.times)
                        Tstress = transport(HEART_TRANSPORT, F.stress; times = HEART_PHANTOM.times)

                        # 4. Phantom: one folder per scenario, every file on one grid
                        scans = [("rest", F.rest, Trest), ("stress", F.stress, Tstress)]
                        R = phantom(HEART_PHANTOM, scans, classify, Vasculature.VoxelGrid(origin, 0.2, dims); F.pairs, F.terr, nvox,
                            out = "heart_healthy", roi, roi_names, roi_groups = ["LV" => [1, 2, 3], "RV" => [4, 5, 6]],
                            tracking = ["descending aorta" => (x, y, z)])"""
                    )
                ),
                P(
                    :class => prose,
                    "Each spec is a plain struct. A variant of an organ is the same call with a keyword changed, ",
                    "for example ", Code(:class => inline, "FlowSpec(HEART_FLOW; p_in = 70)"), "."
                ),

                H2(:id => "load", :class => h2, "Load the Phantom"),
                P(
                    :class => prose,
                    "The scenario folder loads at any time point as material maps for a CT simulator such as ",
                    A(
                        :href => "https://github.com/MolloiLab/BasisSimulator.jl", :target => "_blank",
                        :class => "tv-link", "BasisSimulator.jl"
                    ),
                    "."
                ),
                Pre(
                    :class => code_block,
                    Code(
                        :class => "language-julia", """import BasisSimulator

                        ph = Vasculature.DynamicPhantom.load("heart_healthy"; rois = true) # the bundle
                        P, m = Vasculature.DynamicPhantom.bs_phantom(ph, BasisSimulator; state = "rest", t = 25) # material map at t = 25 s"""
                    )
                ),

                H2(:id => "whats-next", :class => h2, "What's Next"),
                P(
                    :class => prose,
                    "Worked examples are on the ",
                    A(
                        :href => "$(BASE)/examples/",
                        :class => "tv-link", "Examples"
                    ),
                    " page. Every exported function and spec is listed on the ",
                    A(
                        :href => "$(BASE)/api/",
                        :class => "tv-link", "API"
                    ),
                    " page."
                )
            )
        )
    end
end
