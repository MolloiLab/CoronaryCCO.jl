let BASE = get(ENV, "VASCULATURE_BASE", "")
    () -> begin
        prose = "tv-body text-[17px]"
        h2 = "tv-title text-3xl pt-10"
        h3 = "tv-title text-xl pt-4"

        p(xs...) = P(:class => prose, xs...)
        c(s) = Code(:class => "tv-code", s)
        b(s) = Strong(:class => "text-[var(--text)]", s)
        link(href, s) = A(:href => href, :class => "tv-link", s)
        julia(s) = Pre(:class => "tv-pre", Code(:class => "language-julia", s))
        shell(s) = Pre(:class => "tv-pre", Code(:class => "language-bash", s))
        note(title, xs...) = Details(:class => "tv-note", Summary(title), Div(:class => "tv-note-body space-y-3 tv-body text-[15px]", xs...))
        # a rough edge: what you see, why, and what to do
        edge(symptom, cause, fix) = Details(
            :class => "tv-note tv-edge", Summary(symptom),
            Div(
                :class => "tv-note-body space-y-2 tv-body text-[15px]",
                P(b("Why. "), cause...), P(b("What to do. "), fix...),
            ),
        )
        table(head, rows) = Div(
            :class => "tv-table-wrap",
            Table(
                :class => "tv-table",
                Thead(Tr([Th(h) for h in head]...)),
                Tbody([Tr([Td(x) for x in row]...) for row in rows]...),
            ),
        )
        steps(items...) = Ol(:class => "tv-steps space-y-3 tv-body text-[17px]", [Li(x...) for x in items]...)
        bullets(items...) = Ul(:class => "list-disc pl-6 space-y-2 tv-body text-[17px]", [Li(x...) for x in items]...)
        figure(src, alt, caption) = Figure(:class => "tv-figure", Img(:src => "$(BASE)/assets/$src", :alt => alt), Figcaption(caption))

        sections = [
            ("why", "Why Grow a Tree?"),
            ("one-path", "One Code Path"),
            ("install", "Install"),
            ("workflow", "How We Work"),
            ("flow", "Where the Flow Comes From"),
            ("heart", "Walkthrough: the Heart"),
            ("organs", "Three Organs, One Recipe"),
            ("new-organ", "A New Organ"),
            ("new-disease", "A New Disease"),
            ("rough-edges", "Rough Edges"),
            ("rules", "The Rules We Keep"),
            ("next", "Where to Go Next"),
            ("agents", "For Your Agent"),
        ]

        PageWithTOC(
            sections, Div(
                :class => "max-w-3xl space-y-6 pb-24",
                H1(:class => "tv-display text-5xl sm:text-6xl pb-2", "Getting Started"),
                P(
                    :class => "tv-lede",
                    "A CT scan shows an organ's large vessels and nothing smaller. Vasculature.jl grows the rest of the ",
                    "tree, down to the terminal arterioles, then lets physics decide how blood moves through it. It sends a ",
                    "contrast injection through the whole body and into that tree, and writes a CT phantom in which the ",
                    "true perfusion of every voxel is known."
                ),
                p(
                    "This page takes you from installing the package to your first phantom. It shows how the heart, the ",
                    "lungs and the leg were built, and where to start for a new organ or a new disease. It also lists what ",
                    "went wrong for us along the way, so it does not go wrong for you."
                ),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "why", :class => h2, "Why Grow a Tree?"),
                p(
                    "CT perfusion measures blood flow from how fast contrast fills the tissue. To test a method like that, ",
                    "you need tissue whose true flow you know, and patients never tell you. Digital phantoms such as XCAT ",
                    "draw an organ's anatomy beautifully, but they stop at the vessels a scanner resolves: in the heart, ",
                    "the coronary arteries and veins on its surface."
                ),
                p(
                    "So we grow the rest. From XCAT's coronary arteries the heart grows 957,045 terminals, one for every ",
                    "small region of myocardium a terminal arteriole feeds. Flow then follows from first principles. ",
                    "Nothing in the pipeline is told what the perfusion should be: in the healthy heart it comes out at ",
                    "0.86 mL/min/g at rest and 3.85 under adenosine, a flow reserve of 4.5. A CT simulator such as ",
                    link("https://github.com/MolloiLab/BasisSimulator.jl", "BasisSimulator.jl"),
                    " scans the phantom, and every measurement can be checked against the truth, voxel by voxel."
                ),
                figure("step1_grow.png", "The imaged coronary arteries and the tree grown from them", "The coronary arteries XCAT draws (left), and the tree grown from them (right)."),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "one-path", :class => h2, "One Code Path for Every Organ"),
                p("Four functions carry an organ from its imaged vessels to a phantom:"),
                table(
                    ["Stage", "Call", "You give", "You get"],
                    [
                        ["Grow", c("grow(spec, domain, inlets)"), "the tissue as a voxel mask, the imaged vessels as trees, what to avoid", "the grown tree"],
                        ["Flow", c("flow(spec, arteries, veins, tissue)"), "both trees and the tissue's voxel centres", "every vessel's flow at rest and at maximal dilation, every voxel's true perfusion"],
                        ["Transport", c("transport(spec, circuit, injection)"), "a solved circuit and an injection", "the iodine in every vessel and capillary bed, every second"],
                        ["Phantom", c("phantom(spec, scans, classify, grid; …)"), "what each voxel is, and the compartments to report", "one folder: labels, compartments, curves, ground truth"],
                    ],
                ),
                p(
                    "Each stage takes a spec: ", c("GrowthSpec"), ", ", c("FlowSpec"), ", ", c("TransportSpec"), ", ",
                    c("PhantomSpec"), ". A spec is a plain Julia struct, and every value in it is a measurement with a source. ",
                    "The heart is ", c("HEART_GROWTH"), ", ", c("HEART_FLOW"), ", ", c("HEART_TRANSPORT"), " and ",
                    c("HEART_PHANTOM"), "; the lungs are the same four structs with the lungs' numbers. A variant is one ",
                    "keyword: ", c("FlowSpec(HEART_FLOW; p_in = 70)"), "."
                ),
                p(
                    "Nothing in the package asks which organ it is running. A brain will not be a new code path. It will be ",
                    "four new specs, and a notebook that prepares the brain's geometry. ",
                    b("This is the most important design decision in the package, and the one to protect when you add to it.")
                ),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "install", :class => h2, "Install"),
                p("Julia 1.12 or later. Growth runs on the CPU or on any GPU."),
                julia(
                    """using Pkg
                    Pkg.add(url = "https://github.com/MolloiLab/Vasculature.jl")
                    Pkg.add("Metal") # Apple silicon; CUDA for NVIDIA, AMDGPU for AMD, or nothing for the CPU"""
                ),
                p(
                    "The XCAT phantom is licensed. Its files, and everything made from them (meshes, label volumes, grown ",
                    "trees, phantoms, scans), stay on lab storage. They never go into a public repository, an issue or a chat."
                ),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "workflow", :class => h2, "How We Work: Notebooks in SpaceStation"),
                p(
                    "Every organ is built in a Pluto notebook in ", c("docs/notebooks/"), ": ", c("00_cardiac_pipeline"), ", ",
                    c("01_pulmonary_pipeline"), ", ", c("02_limb_pipeline"), ". A notebook is a plain Julia file. It prepares the ",
                    "organ's geometry, calls the four stages, checks the results and draws them. The same file is the ",
                    "documentation: the ", link("$(BASE)/examples/", "Examples"), " pages are these notebooks."
                ),
                p("We run them in SpaceStation, a workspace for Pluto notebooks that people and agents share:"),
                shell("spacestation --agents-md docs/notebooks"),
                p(
                    "It starts in ", b("lazy mode"), ". Nothing runs until you ask, and when a file changes, the cells it ",
                    "touches turn stale instead of running. Outputs are cached beside the notebook and survive restarts. That ",
                    "matters here, because a growth can take minutes and a phantom ten: you run what changed, and nothing else."
                ),
                p(
                    "To work with an agent, start Claude Code in a terminal inside SpaceStation. It edits the notebook as a ",
                    "file, you watch the changed cells go stale, and it runs only those with ",
                    c("pluto-collab run <notebook>.jl --stale"), ". It never reruns the whole notebook, and ",
                    c("pluto-collab restart"), " is only for a worker that died."
                ),
                p(
                    "Every heavy step in a notebook goes through ", c("cached(name) do … end"), ": the first run computes ",
                    "and saves, and every run after loads in seconds. A notebook includes like any Julia file, so a script ",
                    "on a cluster can ", c("include"), " it and reuse every function and cache in it."
                ),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "flow", :class => h2, "Where the Flow Comes From"),
                p(
                    "The rule that shapes everything else is that ", b("flow is never prescribed"), ". Growth knows nothing ",
                    "about flow targets. The flow model knows two pressures, the trees' geometry and the physiology of a vessel ",
                    "wall, and the ground truth is whatever comes out. Here is how the heart gets its flow."
                ),
                steps(
                    (
                        b("Maximal dilation is calibrated once. "),
                        "Under adenosine every arteriole is fully relaxed, so resistance alone sets the flow. Grown arterioles ",
                        "are longer and thinner than real ones (growth runs each terminal to the middle of its territory, where ",
                        "real tissue has a capillary network), so their raw resistance is too high. One factor, κ, scales the ",
                        "microvessels so that the whole bed's resistance matches healthy people measured by PET (",
                        c("r_total_dilated"), "). Their resistances relative to each other, and so the spread of flow from ",
                        "voxel to voxel, still come from the tree.",
                    ),
                    (
                        b("Rest follows the oxygen balance. "),
                        "At rest the heart sets its flow to its oxygen need. Every arteriole below 150 µm is a wall in force ",
                        "balance. Its tone follows wall tension (the myogenic response) and the oxygen state of the beds it ",
                        "feeds (ATP released by desaturating red cells, conducted upstream). One tone constant is set so the ",
                        "blood leaving the heart has the measured oxygen extraction, 0.73. After that, each arteriole finds its ",
                        "own diameter.",
                    ),
                    (
                        b("Everything else is a prediction. "),
                        "The distribution of flow, the pressure–flow curve, transit times and the response to a stenosis are not ",
                        "fitted to anything. They are what the physics gives, and the tests check them against experiments.",
                    ),
                ),
                p(
                    "Oxygen enters in two places and only two. The tissue's O₂ demand (", c("mvo2"), ") and extraction (",
                    c("oef"), ") set resting flow through Fick's principle. Fick also sets how many terminals to grow: ",
                    c("flow_ml_min"), ", the organ's flow, over the flow of one terminal. The lungs use neither. They are not ",
                    "fed by the blood they carry, so their tone is passive, their terminal count is their structure (one per ",
                    "acinus), and their resistance is a prediction too."
                ),
                p(
                    "Every value in a spec has its source in ", c("test/reference/"), ", and the tests compare the model with ",
                    "measurements it was never fitted to. When the two disagree, we say so and leave it. The lungs carry ",
                    "3.45 L/min, a quarter below the measured mean, because their grown veins are narrow; we report that ",
                    "rather than tune it away. A phantom tuned to look right cannot test anything."
                ),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "heart", :class => h2, "Walkthrough: the Heart"),
                p(
                    "The whole heart is the ", link("$(BASE)/examples/00_cardiac_pipeline/", c("00_cardiac_pipeline")),
                    " notebook. These are its four calls, with what each one does."
                ),
                H3(:class => h3, "Grow"),
                p(
                    "The tissue is XCAT's ventricular walls on a 1 mm grid. The inlets are XCAT's coronary arteries as trees ",
                    "of centrelines, and the chambers are blockers. Every imaged vessel first gets a short branch into the wall ",
                    "every 2.8 mm (a penetrating stub), so growth starts along its whole length and not only at its few imaged ends."
                ),
                julia(
                    """using Vasculature, Metal

                    spec = GrowthSpec(HEART_GROWTH; device = MtlArray)
                    tree = grow(spec, "heart.vol", ["lca.xml", "rca.xml"]; blockers = (chambers, origin_mm, 0.4))
                    veins = grow(GrowthSpec(HEART_VEIN_GROWTH; device = MtlArray), "heart.vol", ["veins.xml"])"""
                ),
                p(
                    c("HEART_GROWTH"), " grows 957,045 terminals in about 40 seconds on a laptop GPU. Growth is HCT-CCO's staged ",
                    "constrained constructive optimization: each new terminal joins the tree where it adds the least blood ",
                    "volume, so the tree fills the tissue the way real trees do."
                ),
                H3(:class => h3, "Flow"),
                julia(
                    """F = flow(HEART_FLOW, vtree(tree), vtree(veins), tissue) # tissue: 3 × n voxel centres (mm)
                    F.rest.resist.mbf, F.stress.resist.mbf # mL/min/g at rest and at maximal dilation"""
                ),
                H3(:class => h3, "Transport"),
                p(
                    "An injection into the right arm travels through a whole-body circulation (the venae cavae, the right ",
                    "heart, the lungs, the left heart and every organ bed), into the grown coronary tree, through every ",
                    "capillary bed and vein, and around again."
                ),
                julia(
                    """injection = Body.Injection(; phases = [(35, 5, 350), (35, 5, 105)]) # mL, mL/s, mg I/mL: contrast, then a 30:70 chaser
                    T = transport(HEART_TRANSPORT, F.rest, injection; times = HEART_PHANTOM.times)"""
                ),
                H3(:class => h3, "Phantom"),
                p(
                    "Every 0.2 mm voxel of XCAT's chest is classified: a tissue voxel, a vessel segment, a body pool (the aorta, ",
                    "a ventricle) or a static material. Each voxel also gets a compartment, the region ground truth is reported for."
                ),
                julia(
                    """phantom(HEART_PHANTOM, [("rest", F.rest, Trest), ("stress", F.stress, Tstress)], classify, grid;
                        F.pairs, F.terr, nvox, out = "healthy", roi, roi_names, roi_groups,
                        tracking = ["descending aorta" => (x, y, z)])

                    ph = DynamicPhantom.load("healthy"; rois = true)
                    P, materials = DynamicPhantom.bs_phantom(ph, BasisSimulator; state = "rest", t = 25)"""
                ),
                p(
                    "The folder holds every voxel's label and compartment, each label's iodine every second, the ground-truth ",
                    "perfusion map, and an iodine ledger that must balance: what was delivered equals what is held plus what ",
                    "left through the veins. The build fails if it does not."
                ),
                note(
                    "What a compartment's ground truth means",
                    p(
                        "A compartment is a region of the organ with everything in it: tissue, and every artery, vein and lumen ",
                        "that lies in it. Its ground truth is the flow its capillary beds receive divided by the mass of everything ",
                        "in it, which is what microspheres measure. Groups add compartments together. The heart reports ",
                        "\"LV myocardium\" and \"LV myocardium + epicardial vessels\": a perfusion measurement over the ",
                        "second also sees the iodine still on its way in and out through the surface vessels."
                    ),
                ),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "organs", :class => h2, "Three Organs, One Recipe"),
                p(
                    "The heart, the lungs and the leg use the same four calls. What changes is the geometry each notebook ",
                    "prepares and the numbers in each spec."
                ),
                table(
                    ["", "Heart", "Lungs", "Leg"],
                    [
                        ["Imaged inlets", "coronary arteries and veins", "main pulmonary artery, pulmonary veins", "femoral artery and vein"],
                        ["Terminal count", "Fick: flow over one terminal's flow (957,045)", "structure: one per acinus (29,912)", "Murray: femoral radius and 10 µm terminals (58.7 million)"],
                        ["Bifurcation exponent γ", "2.33", "2.33", "3 (Murray)"],
                        ["Blockers", "the chambers", "airways, fissures, the heart, the aorta", "bone"],
                        ["Tone", "O₂-regulated", "passive", "O₂-regulated"],
                        ["Cardiac output", "the body's", "the lungs' own flow", "the body's"],
                        ["Special step", "penetrating stubs", "lobes from the airway tree", "growth in two passes"],
                        ["Ideal compartment", "myocardium + epicardial vessels", "lungs + hilar vessels", "limb + conduit vessels"],
                    ],
                ),
                p(
                    "The leg shows what happens when an organ is large. Murray's law says a tree's terminal radius is ",
                    c("r_root · N^(−1/3)"), ", so reaching 10 µm from a 3.89 mm femoral artery takes 58.7 million terminals, ",
                    "too many for one growth. So the leg grows twice: first the whole leg down to 10 mL patches, then each ",
                    "patch down to 10 µm on its own. Every radius in both passes follows Murray's law from the terminals it ",
                    "feeds, so the result is one tree. The patches share no vessels, so they run in parallel."
                ),
                figure("step6_limb_two_pass.png", "The leg's two-pass growth", "The leg, a region of ten patches stitched onto the first pass, and 8 mm of it down to 10 µm (02_limb_pipeline)."),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "new-organ", :class => h2, "Building a New Organ"),
                p("Say you want a brain. Here is the order to work in, and where each piece lives."),
                steps(
                    (
                        b("Start the notebook. "), "Copy the closest one (the heart, for an O₂-regulated organ) to ",
                        c("docs/notebooks/03_brain_pipeline.jl"), ", keep its sections, and open the folder in SpaceStation.",
                    ),
                    (
                        b("The tissue. "), "A voxel mask of the tissue to perfuse, written with ", c("write_vol"), " (or ",
                        c("mesh_domain"), " from a surface). For a brain, grey and white matter. A ", c("FlowSpec"),
                        " describes one tissue, so decide early whether they are one or two.",
                    ),
                    (
                        b("The inlets. "), "The vessels the scan shows, as trees: ", c("XCAT.centerline"),
                        " turns a surface into a centreline, ", c("XCAT.base_forest"), " or ", c("XCAT.tube_forest"),
                        " joins them, and ", c("XCAT.write_base_gxl"), " writes them for ", c("grow"),
                        ". For a brain, both internal carotids and both vertebral arteries. Add penetrating stubs (",
                        c("XCAT.penetrating_stubs"), "), and let the imaged lumens be crossed without taking terminals (",
                        c("XCAT.vessel_tube"), " with ", c("NOTERM"), ").",
                    ),
                    (b("The blockers. "), "Anything a vessel must not cross. For a brain, the ventricles and the skull."),
                    (
                        b("The count. "), "Fick for an organ its blood feeds (", c("flow_ml_min"),
                        ": an older adult's 510 mL/min gives 3.6 million terminals), a structural count if the organ's ",
                        "units are known, or Murray from a target terminal size. If it is too large for one growth, grow in ",
                        "two passes, as the leg does.",
                    ),
                    (
                        b("The growth spec. "), c("GrowthSpec(HEART_GROWTH; tissue_ml, root_caps, gamma, flow_ml_min)"),
                        ". Choose γ for the organ, then check it: the imaged vessels should come out near their imaged radii, ",
                        "and no capillary bed should end up with its venous pressure above its arterial one.",
                    ),
                    (
                        b("The flow spec. "), c("FlowSpec(HEART_FLOW; …)"), " with the organ's own values: pressures, blood ",
                        "volume and its split, capillary size, O₂ demand and extraction, resistance at maximal dilation, ",
                        "permeability. Write each source into ", c("test/reference/"), ". Then solve and check: no reversed ",
                        "beds, resting flow and venous saturation within the measured range, tone not saturated.",
                    ),
                    (
                        b("Transport. "), "The body model has organ slots for the heart and the lungs. A brain, a liver or a ",
                        "leg first needs its own slot, with its share of the cardiac output taken from the grown tree's flow, ",
                        "as ", c("coronary_share"), " and ", c("organ_flow"), " do. That is package work: open an issue first.",
                    ),
                    (
                        b("The phantom. "), "A classifier that says what every voxel is, compartments (each also grouped with ",
                        "the conduit vessels that feed and drain it), and a bolus-tracking point inside an imaged artery.",
                    ),
                    (
                        b("Into the package. "), "When the notebook works, move its specs into ", c("src/organs.jl"), " and ",
                        c("src/perfusion/organs.jl"), ", move any helper the next organ would need into ", c("src/"),
                        ", add tests, and add the notebook to the docs.",
                    ),
                ),
                note(
                    "The liver needs one thing the package does not have yet",
                    p(
                        "Two blood supplies. The hepatic artery and the portal vein both feed the sinusoids, so the liver needs two ",
                        "inlet trees, each with its own pressure, in one circuit. Plan that change with the package's maintainers ",
                        "before starting."
                    ),
                ),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "new-disease", :class => h2, "Building a New Disease"),
                p(
                    "A disease is a change to a healthy organ, solved again with the healthy organ's tone constant so its ",
                    "arterioles respond on their own. The heart's 90% LAD stenosis is the example:"
                ),
                julia(
                    """lesion = [(segment, area_fraction, length_mm) for …] # the first 10 mm of the LAD
                    stress = Flow.circuit(arteries, F.terr; state = :stress, veins, vleaf = F.vleaf, phys = HEART_FLOW,
                        reference = F.stress, kappa = F.stress.resist.kappa, stenosis = lesion)
                    rest = Flow.circuit(arteries, F.terr; state = :rest, veins, vleaf = F.vleaf, phys = HEART_FLOW,
                        reference = F.rest, kappa = F.stress.resist.kappa, stenosis = lesion, demand = F.demand)"""
                ),
                p(
                    "The lesion is graded on XCAT's imaged lumen. Its pressure loss follows the viscous and separation terms of ",
                    "Young and Tsai, and the arterioles downstream dilate until they run out of reserve. The phantom then paints ",
                    "plaque (fibrous tissue, a lipid core, calcium) where the lesion narrows the lumen, and the scenario is saved ",
                    "as its own folder."
                ),
                p(
                    "The lungs' occlusion is the same call with a tiny area fraction (10⁻³) on the artery that feeds most of one ",
                    "lobe. More stenoses, serial lesions, microvascular disease or a lower pressure are each a geometry change or ",
                    "a spec field on the same calls. Keep the healthy run as the reference, and give each scenario its own folder."
                ),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "rough-edges", :class => h2, "Rough Edges, and How We Got Past Them"),
                p(
                    "Each of these cost us hours or days. They are written as what you will see, why it happens, and what to do. ",
                    "Open the ones that match what you are looking at."
                ),
                H3(:class => h3, "Geometry and inputs"),
                edge(
                    "The tree grows only from the imaged vessels' ends",
                    ("Staged growth attaches new vessels near the inlet tree's ends. An imaged artery drawn without its side branches gives growth almost nowhere to start.",),
                    ("Add penetrating stubs (", c("XCAT.penetrating_stubs"), "): a short branch into the tissue about every 2.8 mm × the organ's size scale, wherever tissue is in reach."),
                ),
                edge(
                    "New branches collide with the next imaged piece",
                    ("A branch leaves its host on the host's centreline, and the collision test exempts only the host. Imaged nodes closer together than the host's radius make every branch hit its neighbour.",),
                    ("Resample imaged centrelines every 2.8 mm × the organ's size scale before joining them into trees.",),
                ),
                edge(
                    "Imaged vessels run outside the tissue",
                    ("XCAT draws the epicardial arteries, the hilar vessels and the femoral artery partly outside the tissue mask, so growth cannot leave them.",),
                    ("Add their lumens (and short bridges to the tissue) to the domain as ", c("NOTERM"), " voxels: crossable, but no terminal ends there (", c("XCAT.vessel_tube"), ", ", c("XCAT.vessel_bridges"), ")."),
                ),
                edge(
                    "XCAT's anatomy is outside the measured range",
                    ("XCAT vmale_50's aorta is about three standard deviations narrow: 21–22 mm ascending against 33 ± 4 mm in adult men. Its blood arrived too early. Its main pulmonary artery is large, and its femoral vein narrow.",),
                    ("Check every vessel you use against measurements. For positional great vessels, take the timing from the body model's volumes (", c("XCAT.scale_tree"), "), and write down what you kept as drawn."),
                ),
                edge(
                    "A vessel is labelled as part of a chamber",
                    ("XCAT's label volume paints the main pulmonary artery as the right ventricle's blood.",),
                    ("Use the vessel's NURBS surface, not its label, to tell the two apart.",),
                ),
                edge(
                    "A label means something else than you think",
                    ("XCAT reuses label ids in distant organs, its airway label 7 is the wall (cartilage) and 11 the air inside it, and its surface meshes sit at a fixed shift from its NURBS and label frame.",),
                    ("Restrict every label to the region you expect it in, look each id up before using it, and apply the mesh shift once, where the volume is built.",),
                ),
                edge(
                    "Cutting the imaged arteries short so growth can define them",
                    ("Below about 60 mm of imaged artery, growth stalls; the large vessels must come from the image.",),
                    ("Keep the imaged trees whole. Where tissue lies far from every imaged vessel (15.5% of one lung lay over 20 mm away), add supply branches (", c("XCAT.fill_supply_gaps"), ")."),
                ),
                edge(
                    "A new organ's root radius is the heart's",
                    ("With ", c("root_caps = nothing"), ", growth falls back to 1.588 mm, a coronary value."),
                    ("Always set ", c("root_caps"), " to the imaged inlet radius."),
                ),
                H3(:class => h3, "Growth"),
                edge(
                    "Capillary beds reverse, or the conduits take all the pressure",
                    ("The bifurcation exponent γ is organ-specific. Grown with the coronary 2.33, the leg's conduits came out at 0.54–0.76 of their imaged radii and 39% of its beds reversed.",),
                    ("Choose γ per organ and check the imaged vessels' grown radii against their imaged ones. The leg uses 3.",),
                ),
                edge(
                    "Terminals are far larger than real terminal arterioles",
                    ("A tree's terminal radius is r_root · N^(−1/γ). Too few terminals means terminals that are too big. Lowering γ shrinks them, but starves them of pressure: about 800 mmHg lost in a 10 mL block.",),
                    ("Grow more terminals, never a lower γ. If that is too many for one growth, grow in two passes (", c("02_limb_pipeline"), ")."),
                ),
                edge(
                    "Second-pass patches end in terminals that are too thin",
                    ("Growth balances each terminal's path resistance, so first-pass terminals that carry the same flow scatter ±35% in radius. A patch grown from a thin one ended at 6 µm and lost 128 mmHg.",),
                    ("Take every radius from Murray's law and the final terminal counts below it, r = 10 µm · n^(1/3), in both passes.",),
                ),
                edge(
                    "Seams between patches",
                    ("Growth keeps terminals a clearance away from the domain's edge, so each patch's border is left empty.",),
                    ("Grow each patch with a 3 mm margin of its neighbours' tissue as ", c("NOTERM"), "."),
                ),
                edge(
                    "Heart-sized constants misbehave in a bigger organ",
                    ("Distances inside growth (how far a new terminal may connect, the stub reach, the search grid) were first sized for a 148 mL heart.",),
                    ("Scale them with the organ's size (", c("size_scale(spec)"), ") and check each one is a numerical guard, not physiology: changing it should not change the tree."),
                ),
                H3(:class => h3, "Flow"),
                edge(
                    "NaN pressures, or reversed beds, on a new organ",
                    ("The heart's ", c("FlowSpec"), " asks a leg for the heart's flow per gram: 47 L/min for one leg."),
                    ("Give every organ its own ", c("FlowSpec"), ", every value measured in that organ and cited."),
                ),
                edge(
                    "Flow comes out below the measured mean",
                    ("The flow is a prediction from the pressures and the trees, and a tree can be imperfect: the lungs' narrow grown veins give 3.45 L/min against 4.8.",),
                    ("Report it, with the reason. Do not tune growth or the spec to hit the number.",),
                ),
                edge(
                    "Perfusion at stress is measured far below the truth",
                    ("This was the physics, not the measurement. Contrast permeability was a constant taken from a Patlak influx constant (a lower bound), and adenosine's rise in blood volume was missing.",),
                    ("Permeability now rises with flow (Renkin–Crone, PS = 0.22 + 0.35 · MBF) and the dilated blood volume is measured (", c("fb_dilated"), "). Fix a disagreement in the physics, never in the measurement method."),
                ),
                edge(
                    "A total occlusion gives NaN",
                    ("A segment with zero area has infinite resistance.",),
                    ("Use a tiny area fraction instead (the lungs' occlusion uses 10⁻³).",),
                ),
                edge(
                    "Imaged arteries constrict at rest",
                    ("Tone applies to every artery whose radius is below the active cut-off, imaged or not, and a grown radius can be small even where the image is large.",),
                    ("That is intended: an imaged vessel is taken at the size the tree gives it. Check the grown radii of imaged vessels against the image (see the γ edge above).",),
                ),
                H3(:class => h3, "Transport and timing"),
                edge(
                    "The body and the organ carry different flows",
                    ("The whole-body model and the grown organ are solved separately.",),
                    ("The body's share for the organ comes from the grown tree: ", c("coronary_share"), " for the heart, ", c("organ_flow"), " for the lungs (which carry the whole cardiac output). A new organ needs the same."),
                ),
                edge(
                    "The arterial input and the tissue are out of step",
                    ("An arterial curve measured upstream of a compartment's inlet arrives early; one measured in a large conduit carries iodine the tissue has not yet seen.",),
                    ("Time each great vessel from its own position along the body. Report the ideal compartment with its conduit vessels included (epicardial, hilar), so iodine in transit is counted.",),
                ),
                edge(
                    "The injection is not the one you meant",
                    ("", c("Body.Injection()"), " with no arguments is 50 mL at 5 mL/s and a 40 mL saline flush, not the lab's protocol."),
                    ("Pass the injection explicitly, the same in every organ: 0.5 mL/kg of contrast at 5 mL/s, then 0.5 mL/kg of a 30:70 contrast:saline chaser.",),
                ),
                edge(
                    "The arterial peak is lower at stress",
                    ("Adenosine raises cardiac output 1.71×, so the same injection is diluted into more blood.",),
                    ("Nothing: it is real indicator dilution.",),
                ),
                H3(:class => h3, "Phantom and painting"),
                edge(
                    "The bolus-tracking curve differs between scenarios, or reads low",
                    ("Vessels are grouped into classes by their curves, and an imaged vessel's iodine is spread over XCAT's lumen, which can be wider than the grown radius. The pulmonary trunk's curve read 19% below its blood.",),
                    ("A tracking point's segment gets a class of its own, painted at its true blood concentration (", c("phantom"), " does this). Keep a CT arterial ROI inside that segment."),
                ),
                edge(
                    "The iodine ledger fails",
                    ("Painted iodine must equal the iodine the transport put in the organ, within 1% of the peak delivered.",),
                    ("A failing ledger is a classifier or compartment bug. Check which voxels a segment paints before anything else.",),
                ),
                edge(
                    "Too many materials for the CT simulator, or too much memory to load",
                    ("BasisSimulator takes at most 64 materials, and a full-chest bundle is 4.5 billion voxels.",),
                    ("Merge labels into materials per time point with ", c("DynamicPhantom.materials_at"), " (every compartment keeps its iodine), and load a slab with ", c("DynamicPhantom.load(dir; zrange)"), "."),
                ),
                H3(:class => h3, "Notebooks and workflow"),
                edge(
                    "A loop in a notebook takes minutes",
                    ("Loops over a notebook's global variables are not compiled for their types.",),
                    ("Put the loop in a function and pass the data in (a function barrier).",),
                ),
                edge(
                    "A figure with thousands of vessels takes minutes",
                    ("One ", c("lines!"), " per vessel means one plot object per vessel."),
                    ("Draw them all with one ", c("linesegments!"), "."),
                ),
                edge(
                    "A cache fails to load after a package change",
                    ("A cached result holds the old struct, for example a ", c("FlowSpec"), " before a field was added."),
                    ("Delete that cache file, and let the notebook rebuild it.",),
                ),
                edge(
                    "An edited cell loses its last character",
                    ("Pluto's file format expects a blank line between cells, and drops the last character of a cell that has none.",),
                    ("Keep a blank line before every ", c("# ╔═╡"), " marker, and check ", c("pluto-collab status"), " after an edit."),
                ),
                edge(
                    "A notebook cell will not parse, or a name is ambiguous",
                    ("A docstring cannot sit before a ", c("begin"), " or ", c("let"), " block, a one-line ", c("do … end"), " does not parse, ", c("Slider"), " is exported by both PlutoUI and CairoMakie, and ", c("map"), " does not take ", c("pairs"), " of a NamedTuple."),
                    ("Use a comment instead of the docstring, write the block on several lines, import names explicitly (", c("using PlutoUI: Slider"), "), and use a comprehension."),
                ),
                edge(
                    "A package change does not show up in the notebook",
                    ("The notebook's worker loaded the package when it started.",),
                    ("Restart the worker after changing the package or its environment. Never kill workers by process name: that also kills the live session's.",),
                ),
                edge(
                    "The docs deploy fails with a stale export",
                    ("A notebook page is fingerprinted by its source. A live SpaceStation session re-saves the file, so the committed source can differ from the one that was rendered.",),
                    ("Render (", c("docs/extract_all.jl"), "), check the notebook file did not change, verify (", c("docs/verify_notebook_exports.jl"), "), and commit both together. Any change in ", c("src/"), " re-renders every notebook page."),
                ),
                edge(
                    "A new notebook breaks the docs build",
                    ("The docs site lists notebooks it knows, and the verifier rejects machine paths (", c("/Users"), ", ", c("/Volumes"), ") and errors in a page's outputs. A notebook with sliders would have every slider setting precomputed."),
                    ("Add the notebook to ", c("NOTEBOOK_META"), " and the page titles, list its sliders in ", c("FORCE_FALLBACK_BONDS"), " (", c("docs/extract_all.jl"), "), keep local paths out of outputs, and commit its render."),
                ),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "rules", :class => h2, "The Rules We Keep"),
                bullets(
                    (b("One code path. "), "One function per stage, the same for every organ. Organs differ only in their specs and their geometry; nothing in ", c("src/"), " branches on the organ."),
                    (b("Specs are Julia structs. "), "Presets live in code, every value cited, with its source in ", c("test/reference/"), ". No configuration files."),
                    (b("Flow is pressure-driven. "), "Ground truth is what the physics gives. Measurements are comparators, never targets, and a disagreement is reported, not tuned away."),
                    (b("Each organ has its own physiology. "), "Never carry one organ's values, thresholds or exponent into another."),
                    (b("Grow real terminals. "), "Reach terminal-arteriole size with more terminals, never with a lower γ. Growth uses HCT-CCO's staged schedule only."),
                    (b("A disease keeps the healthy tone. "), "Solve it with the healthy circuit as its ", c("reference"), ", so the arterioles respond on their own."),
                    (b("Notebooks orchestrate; the package holds the logic. "), "A helper a second organ needs moves into ", c("src/"), ". Superseded code is deleted, not archived."),
                    (b("Small numerical differences are fine. "), "CPU, Metal and CUDA runs differ in the last digits. Do not chase bit-exact results."),
                    (b("Licensed data stays on lab storage. "), "XCAT, and everything made from it, never goes anywhere public."),
                ),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "next", :class => h2, "Where to Go Next"),
                bullets(
                    ("The ", link("$(BASE)/examples/", "Examples"), ": the heart, the lungs and the leg, each as a notebook you can open in SpaceStation."),
                    ("The ", link("$(BASE)/api/", "API"), ": every exported function and spec, with its docstring."),
                    ("The source: ", c("src/spec.jl"), " (growth), ", c("src/perfusion/flow.jl"), " (flow), ", c("src/perfusion/transport.jl"), " and ", c("src/perfusion/body.jl"), " (transport), ", c("src/phantom.jl"), " (phantom), ", c("src/xcat/"), " (XCAT inputs), ", c("src/organs.jl"), " and ", c("src/perfusion/organs.jl"), " (the presets)."),
                ),

                # ─────────────────────────────────────────────────────────────────
                H2(:id => "agents", :class => h2, "For Your Agent"),
                p(
                    "Paste this into any coding agent before it builds on the package. It holds the goals, the rules and the ",
                    "rough edges above in the form an agent needs."
                ),
                CodeBlock(
                    """You are working on Vasculature.jl (github.com/MolloiLab/Vasculature.jl), a Julia package that grows complete vascular trees from imaged vessels, solves blood flow from pressures, transports contrast through a whole-body model, and writes CT perfusion phantoms with known ground truth. Docs: molloilab.github.io/Vasculature.jl.

                    PIPELINE (one code path for every organ; organs differ only in specs and geometry)
                    - grow(spec::GrowthSpec, domain, inlets; blockers, territories) -> engine; vtree(engine) — src/spec.jl
                    - flow(spec::FlowSpec, arteries, veins, tissue_centres) -> (; rest, stress, terr, pairs, vleaf, demand); for a disease, Flow.circuit(...; reference, kappa, stenosis, demand) — src/perfusion/flow.jl
                    - transport(spec::TransportSpec, circuit, injection; times) — src/perfusion/transport.jl, body.jl
                    - phantom(spec::PhantomSpec, scans, classify, grid; pairs, terr, nvox, out, roi, roi_names, roi_groups, tracking) -> DynamicPhantom bundle — src/phantom.jl
                    - Presets: HEART_*, LUNG_* in src/organs.jl and src/perfusion/organs.jl; sources in test/reference/.
                    - XCAT inputs: src/xcat/ (centerline, base_forest, tube_forest, graft, drain, penetrating_stubs, vessel_tube, vessel_bridges, scale_tree, great_vessel_transit).
                    - Organ notebooks: docs/notebooks/00_cardiac_pipeline.jl, 01_pulmonary_pipeline.jl, 02_limb_pipeline.jl.

                    RULES
                    - Never add organ-specific branches to src/. A new organ = new GrowthSpec/FlowSpec/TransportSpec/PhantomSpec + a notebook that prepares its geometry. Helpers a second organ needs move into src/ with tests.
                    - Specs are Julia structs; every value measured and cited (test/reference/*.toml). No config files.
                    - Flow is never prescribed. Ground truth is what the physics gives. Do not tune growth or specs to hit a measurement; report disagreements with the reason (e.g. lung flow 3.45 vs 4.8 L/min).
                    - Each organ has its own physiology; never reuse another organ's values, thresholds or gamma.
                    - Reach real terminal size (~10 µm) with more terminals (two-pass growth if needed), never by lowering gamma. Growth uses the HCT-CCO staged schedule only.
                    - Diseases: re-solve with the healthy circuit as `reference` (keeps the tone constant); one folder per scenario.
                    - Never change a perfusion measurement method to agree with the phantom; fix the physics.
                    - Delete superseded code; no archives or .bak files. Float noise across CPU/Metal/CUDA is fine.
                    - XCAT data and anything derived from it is licensed: never commit, publish or paste it.

                    WORKFLOW (Pluto notebooks in SpaceStation, lazy mode)
                    - Notebooks are plain .jl files; edit them as files. Keep a blank line before every "# ╔═╡" cell marker (Pluto drops a cell's last character otherwise). Update the "# ╔═╡ Cell order:" block when adding cells; fold table/figure/markdown cells (╟─), keep TableOfContents() shown.
                    - Run only what changed: `pluto-collab run <nb>.jl --stale`; check `pluto-collab status <nb>.jl`. Never rerun a whole notebook; `pluto-collab restart` only for a dead worker.
                    - Wrap heavy steps in cached(name) do ... end; delete a cache after changing a struct it holds.
                    - Put loops over notebook globals in functions; draw many segments with one linesegments!.
                    - Docs: render with docs/extract_all.jl, verify with docs/verify_notebook_exports.jl, and commit the notebook and its export together (a live session can re-save the notebook after rendering).

                    ROUGH EDGES ALREADY SOLVED (do not rediscover)
                    - Growth only from imaged ends -> XCAT.penetrating_stubs (2.8 mm x size scale).
                    - Branches collide with imaged pieces -> resample centrelines every 2.8 mm x size scale.
                    - Imaged vessels outside the tissue -> vessel_tube / vessel_bridges as NOTERM voxels.
                    - XCAT anatomy off-normal (aorta ~3 SD narrow) -> time positional vessels from the body model (scale_tree); check every vessel against measurements.
                    - Beds reverse with the coronary gamma on a large organ -> gamma per organ (leg: 3).
                    - Terminals too large -> more terminals; r_terminal = r_root * N^(-1/gamma).
                    - Two-pass patches too thin / seams -> Murray radii from final counts (r = 10 µm * n^(1/3)); 3 mm NOTERM margin into neighbouring patches.
                    - NaN pressures on a new organ -> the organ's own FlowSpec.
                    - Stress perfusion low -> PS rises with flow (Renkin-Crone, PS = 0.22 + 0.35 MBF), fb_dilated.
                    - Body vs organ flow -> coronary_share (heart), organ_flow (lungs); new organs need a body slot (ORGAN_SLOTS is (:coronary, :lung) today) with the same rule.
                    - AIF timing -> positional great-vessel transit; ideal compartment includes conduit vessels.
                    - Tracking curve wrong -> tracked segment has its own class at blood concentration; AIF ROI inside it.
                    - Ledger must balance within 1%; a failure is a classifier/compartment bug.
                    - XCAT labels: ids reused by distant organs (restrict by region); label 7 = airway wall, 11 = air; the main pulmonary artery is painted as RV blood (use its NURBS surface); meshes are shifted from the NURBS/label frame (apply the shift once).
                    - Never truncate imaged arteries to "let growth define them" (stalls below ~60 mm); fill far tissue with XCAT.fill_supply_gaps.
                    - grow() falls back to root_caps = 1.588 mm (a coronary value): always set the imaged inlet radius.
                    - Size-dependent growth constants (connection cap, stub reach, search grid) scale with size_scale(spec) and must be numerical guards, not physiology.
                    - Total occlusion (area 0) -> NaN; use a tiny area fraction (1e-3).
                    - Body.Injection() default is 50 mL + 40 mL saline, NOT the lab protocol (0.5 mL/kg contrast at 5 mL/s, then 0.5 mL/kg of 30:70 contrast:saline); always pass it explicitly.
                    - Lower arterial peak at stress is real dilution (CO x1.71), not a bug.
                    - BasisSimulator takes <= 64 materials: DynamicPhantom.materials_at; load slabs with zrange.
                    - Notebook traps: no docstring before begin/let; no one-line do...end; Slider is in both PlutoUI and CairoMakie; map does not take pairs(NamedTuple); restart the worker after a package change; never pkill workers by name; CairoMakie needs XDG_CACHE_HOME.
                    - Docs: a new notebook needs NOTEBOOK_META + page title + committed render; sliders go in FORCE_FALLBACK_BONDS; no machine paths (/Users, /Volumes) or errors in outputs; any src/ change re-renders every page.
                    - Public repo: no XCAT data, no XCAT-derived files, and no wording about the lab's private perfusion-analysis code.

                    OPEN WORK
                    - Body slots and flow share for organs other than heart and lungs (brain, liver, leg).
                    - Two blood supplies (liver: hepatic artery + portal vein).
                    - Flow on very large trees (the leg: ~117 M segments) may need a hierarchical solve per patch.""";
                    language = "text", id = "agent-context",
                ),
            )
        )
    end
end
