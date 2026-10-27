let BASE = get(ENV, "VASCULATURE_BASE", ""),
        # One entry per exported name: (id, signature, description). The text
        # is the first sentences of each docstring; the docstrings are the
        # full reference (`?grow` in the REPL).
        API = [
        (
            id = "grow",
            title = "Grow",
            intro = "Grow a vascular tree from the imaged inlets into the organ's domain.",
            entries = [
                (
                    "GrowthSpec(; nterm, flow_ml_min, terminal_flow_ml_min, tissue_ml, reference_ml, wall_clearance, gamma, root_caps, root_fixed, split, symmetry_floor, seed, T, device)",
                    "Everything that differs between organs when growing a vascular tree. Every field is a keyword and every one is required. A run varies an organ with GrowthSpec(spec; kwargs...).",
                ),
                (
                    "HEART_GROWTH, HEART_VEIN_GROWTH, LUNG_GROWTH, LUNG_VEIN_GROWTH",
                    "The organs' GrowthSpec instances.",
                ),
                (
                    "grow(spec::GrowthSpec, domain, inlets; root, blockers, territories, base_groups) -> Engine",
                    "Grow a vascular forest from the imaged inlets (GXL files, one tree each) into the tissue mask domain (a .vol file, 1 mm voxels). blockers is (mask, origin_mm, cell_mm), a fine mask no new vessel may cross.",
                ),
                (
                    "mesh_domain(path, vertices, triangles; margin = 3) -> path",
                    "The growth domain of an organ given as a closed surface mesh: the voxels of a 1 mm grid whose centres lie inside the mesh, written to path as a .vol.",
                ),
                (
                    "terminal_count(spec::GrowthSpec) -> Int",
                    "spec.nterm when given; otherwise the organ's flow over the flow of one terminal.",
                ),
                (
                    "size_scale(spec::GrowthSpec) -> Float32",
                    "The organ's length against the reference, (tissue_ml / reference_ml)^(1/3); exactly 1 for the reference heart.",
                ),
                (
                    "Unified(wrap, copy)",
                    "A device that shares memory with the host (Apple silicon, integrated GPUs). Each batched round hands the growing tree to the device in place instead of copying it.",
                ),
            ],
        ),
        (
            id = "flow",
            title = "Flow",
            intro = "Blood flow in the grown bed from its two pressures alone.",
            entries = [
                (
                    "FlowSpec(; kwargs...)",
                    "Boundary conditions, tissue constants and vessel-wall mechanics of one organ's vascular bed. Every field is a keyword and every one is required. A variant: FlowSpec(HEART_FLOW; p_in = 70).",
                ),
                (
                    "HEART_FLOW, LUNG_FLOW",
                    "The organs' FlowSpec instances, each value commented with its sources.",
                ),
                (
                    "vtree(E) -> VTree",
                    "A grown Engine as a vessel list for flow.",
                ),
                (
                    "flow(spec::FlowSpec, arteries, veins, tissue::AbstractMatrix; to = identity) -> (; rest, stress, terr, pairs, vleaf, terms, vterms, demand)",
                    "Blood flow from the pressures, the vessel sizes and the arterioles' tone. tissue is 3 × n points (mm), one per tissue voxel.",
                ),
                (
                    "couple_veins(arteries, terr, veins, terms, vterms; phys, reach = 3, passes = 24, tol = 5e-4, to = identity)",
                    "Which venule drains each capillary bed, starting from a balanced pairing.",
                ),
            ],
        ),
        (
            id = "transport",
            title = "Transport",
            intro = "A contrast injection through the whole body and the organ's grown tree.",
            entries = [
                (
                    "TransportSpec(; phys, tissue_ml, bed, input, co, beat, f_coronary, tmax, dt, every, coronary)",
                    "Everything that sets where an injected bolus goes. Every field is a keyword and every one is required.",
                ),
                (
                    "HEART_TRANSPORT, LUNG_TRANSPORT",
                    "The organs' TransportSpec instances.",
                ),
                (
                    "transport(spec::TransportSpec, C, injection = Injection(); times = 0:60, ...) -> (; body, input, t_out, outflow, artery, vein, tissue, passes, change)",
                    "The injection through the whole body and the organ's own solved circuit C (from flow). The organ's venous outflow is fed back into the body until it stops changing, so recirculation is included.",
                ),
                (
                    "circulation(spec::TransportSpec)",
                    "The whole-body circulation with the organ's own values in its body segment.",
                ),
                (
                    "coronary_share(C, co) -> Float32",
                    "The coronary bed's share of the cardiac output co (mL/s) when the heart is the grown circuit C.",
                ),
            ],
        ),
        (
            id = "phantom",
            title = "Phantom",
            intro = "The voxelized dynamic phantom of a transported bed.",
            entries = [
                (
                    "PhantomSpec(; times, tissue_material, blood_material, tissue_classes, vessel_classes, vessel_rmin_mm, delay_bin_s, energy_keV, materials, pools)",
                    "Everything that sets an organ's dynamic CT phantom. Every field is a keyword and every one is required.",
                ),
                (
                    "HEART_PHANTOM, LUNG_PHANTOM",
                    "The organs' PhantomSpec instances.",
                ),
                (
                    "phantom(spec, scans, classify, grid; pairs, terr, nvox, out, roi, roi_names, roi_groups, tracking, facts, check = true)",
                    "Writes the phantom of a transported bed into the scenario folder out: labels, ROIs, ground-truth blood flow, per-label and per-ROI iodine curves, the iodine ledger and phantom.toml.",
                ),
                (
                    "VOX_STATIC, VOX_POOL, VOX_TISSUE, VOX_VESSEL",
                    "The voxel kinds classify(p) returns: a static material, a blood pool of the body model, perfused tissue, or a grown vessel.",
                ),
                (
                    "Pool(name, material, segment, what; scale = 1, ml = 0)",
                    "A compartment filled from the whole-body model.",
                ),
                (
                    "pool_curve(P::Pool, r, t; delay = 0, tissue = 0) -> Float32",
                    "Pool P's iodine (mg per mL) at time t (s) in a transport result's body.",
                ),
                (
                    "XCAT_BODY_POOLS, XCAT_STATIC",
                    "The pools and static materials of the XCAT body around the perfused organ.",
                ),
                (
                    "Vasculature.DynamicPhantom.load(dir; volume = true, rois = false, zrange = nothing)",
                    "Loads a scenario folder. zrange loads a large phantom a slab at a time.",
                ),
                (
                    "Vasculature.DynamicPhantom.bs_phantom(ph, BasisSimulator; state, t, max_materials = 64)",
                    "The phantom at time t (s) as a BasisSimulator.jl Phantom.",
                ),
            ],
        ),
        (
            id = "lower-level",
            title = "Lower-Level Tools",
            intro = "The growth engine and tree files underneath grow.",
            entries = [
                (
                    "CCOParams(; kwargs...)",
                    "Run parameters of the growth engine. Defaults match the upstream HctCco command line.",
                ),
                ("CCOModel{T}", "The tree as a struct of arrays: one entry per segment, plus global hydraulic state."),
                ("MaskDomain{T}", "The image-mask domain: the mask, its exact distance map, and the wall clearance."),
                ("Engine{T}", "A growing tree: the model, its spatial grid and its domain. grow returns one."),
                ("grow!(E; parallel = true, progress = 0, growth = :batched, ...)", "Add terminals until nterm."),
                ("nsegments(m)", "The number of segments of a tree."),
                ("set_forbidden!(dom, mask, origin, h)", "A fine forbidden mask no new vessel segment may cross."),
                ("set_territories!(dom, allow)", "Hierarchical territories on the domain grid: which tree groups may take a terminal at each voxel."),
                ("read_vol(path) -> (img, lower)", "Read a DGtal .vol voxel file."),
                ("write_vol(path, mask, lower) -> path", "Write a Bool mask as a DGtal .vol (1 mm voxels)."),
                ("import_gxl(T, path, params) -> CCOModel{T}", "Read a tree (or, with a vector of paths, a forest) from GXL files."),
                ("export_gxl(m, path)", "Write a tree as GXL with flow, resistance and radius per segment."),
                ("dds_orders(m; maxiter = 100, tol = 0.01) -> NamedTuple", "Diameter-defined Strahler ordering (Jiang, Kassab and Fung 1994)."),
                ("tree_qc(E; rmin = 0.25, hops = 2) -> NamedTuple", "Checks a grown tree for overlapping vessels and sharp turns."),
            ],
        ),
    ]

    () -> begin
        card_cls = "tv-tile p-6 sm:p-7 space-y-2"
        sig_cls = "font-mono text-[13px] font-medium break-words"
        prose_cls = "tv-body text-[15px]"
        h2_cls = "tv-title text-3xl pt-6"
        inline = "tv-code"

        sections = [("overview", "Overview"); [(s.id, s.title) for s in API]]

        PageWithTOC(
            sections, Div(
                :class => "max-w-3xl space-y-8",
                H1(
                    :class => "tv-display text-5xl sm:text-6xl",
                    "API Reference"
                ),

                H2(:id => "overview", :class => h2_cls, "Overview"),
                P(
                    :class => prose_cls,
                    "Each stage is one function and one spec struct: ",
                    Code(:class => inline, "grow"), " with ", Code(:class => inline, "GrowthSpec"), ", ",
                    Code(:class => inline, "flow"), " with ", Code(:class => inline, "FlowSpec"), ", ",
                    Code(:class => inline, "transport"), " with ", Code(:class => inline, "TransportSpec"), ", and ",
                    Code(:class => inline, "phantom"), " with ", Code(:class => inline, "PhantomSpec"), ". ",
                    "An organ is one instance of each spec. This page lists the exported names; the docstrings ",
                    "hold the full reference (type ", Code(:class => inline, "?grow"), " in the REPL). The ",
                    A(:href => "$(BASE)/examples/", :class => "tv-link", "Examples"),
                    " show them in a real pipeline."
                ),

                [
                    Div(
                        :class => "space-y-4",
                        H2(:id => s.id, :class => h2_cls, s.title),
                        P(:class => prose_cls, s.intro),
                        [
                            Div(
                                :class => card_cls,
                                Div(:class => sig_cls, sig),
                                P(:class => prose_cls, text)
                            ) for (sig, text) in s.entries
                        ]...
                    ) for s in API
                ]...,

                Div(
                    :class => "pt-6",
                    P(
                        :class => prose_cls,
                        "XCAT inputs (NURBS surfaces, vessel centrelines, organ masks) are read with ",
                        Code(:class => inline, "Vasculature.XCAT"), " (", Code(:class => inline, "read_nrb"), ", ",
                        Code(:class => inline, "centerline"), ", ", Code(:class => inline, "base_forest"), ", ",
                        Code(:class => inline, "penetrating_stubs"), ", ", Code(:class => inline, "fill_supply_gaps"),
                        ")."
                    )
                )
            )
        )
    end
end
