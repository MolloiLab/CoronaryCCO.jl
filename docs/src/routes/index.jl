let BASE = get(ENV, "VASCULATURE_BASE", ""),
        # The four stages, each with its README picture. Defined inside the `let`
        # so Therapy's file-based router doesn't register a route for it.
        STAGES = [
        (
            index = "1",
            title = "Grow",
            image = "step1_grow.png",
            alt = "Imaged input vessels and the tree grown from them",
            caption = "The imaged coronaries (left) and the tree grown from them (right)",
            text = "Start from the vessels a scan shows and grow the rest of the tree inside the organ. " *
                "A million terminals in about 40 seconds on a laptop GPU.",
        ),
        (
            index = "2",
            title = "Flow",
            image = "step2_flow.png",
            alt = "Blood flow computed from pressure",
            caption = "Flow against perfusion pressure, at rest and at maximal dilation, beside human PET",
            text = "Every vessel's flow follows from the inlet and outlet pressures, the vessel sizes and the " *
                "arterioles' tone. No flow is prescribed.",
        ),
        (
            index = "3",
            title = "Transport",
            image = "step3_transport.png",
            alt = "Contrast transport through the body and the grown tree",
            caption = "When the bolus reaches each stage on its way to the heart, and the enhancement it gives",
            text = "A clinical injection travels through the whole body into the grown tree, then through " *
                "every capillary bed and vein, with recirculation.",
        ),
        (
            index = "4",
            title = "Phantom",
            image = "step4_heart_phantom.png",
            alt = "Heart input phantom of XCAT vmale_50",
            caption = "The XCAT chest at 0.2 mm, axial and coronal, coloured by voxel class",
            text = "One folder per scenario on one grid: labels, ROIs, ground-truth blood flow and iodine " *
                "every second, ready for a CT simulator.",
        ),
    ],
        MOVIES = [
        ("heart_flow.webp", "Heart", "contrast through the grown coronary tree"),
        ("lung_flow.webp", "Lungs", "contrast through the grown pulmonary trees"),
    ]

    () -> begin
        Div(
            # ─── Hero: calm type over a faint ambient glow ─────────────────────
            Div(
                :class => "relative",
                Div(:class => "tv-glow"),
                Div(
                    :class => "relative z-10 max-w-3xl space-y-7 pt-6 sm:pt-12",
                    P(:class => "tv-eyebrow", "Grow · Flow · Transport · Phantom"),
                    H1(
                        :class => "tv-display text-5xl sm:text-6xl lg:text-7xl",
                        "Perfusion Phantoms from Grown Vascular Trees"
                    ),
                    P(
                        :class => "tv-body text-lg sm:text-xl max-w-2xl",
                        "Grow a complete vascular tree inside any organ from the few vessels a scan shows, compute its " *
                            "blood flow from the pressures alone, send a contrast injection through it, and export a " *
                            "voxelized perfusion phantom for CT simulation. One algorithm for every organ; an organ " *
                            "differs only in the numbers you give it."
                    ),
                    Div(
                        :class => "flex flex-wrap gap-3 pt-1",
                        A(:href => "$(BASE)/getting-started/", :class => "tv-pill tv-pill-on", "Get Started"),
                        A(:href => "$(BASE)/examples/", :class => "tv-pill", "Examples"),
                        A(
                            :href => "https://github.com/MolloiLab/Vasculature.jl", :target => "_blank",
                            :class => "tv-pill", "GitHub"
                        )
                    )
                )
            ),

            # ─── Heart and lung: two wide tiles ────────────────────────────────
            Div(
                :class => "relative z-10 mt-16 sm:mt-20 grid grid-cols-1 md:grid-cols-2 gap-5 sm:gap-6",
                [
                    Figure(
                        :class => "tv-tile tv-enter", :style => "animation-delay:$(80 * (k - 1))ms",
                        Div(
                            :class => "tv-tile-media",
                            Img(:src => "$(BASE)/assets/$(file)", :alt => "$(label): $(alt)", :class => "w-full h-auto block")
                        ),
                        Figcaption(:class => "tv-tile-label", Strong(label), " · $(alt)")
                    ) for (k, (file, label, alt)) in enumerate(MOVIES)
                ]...
            ),

            # ─── The four stages: big rounded tiles ────────────────────────────
            Div(
                :class => "mt-28 sm:mt-36",
                Div(
                    :class => "max-w-2xl mb-10 space-y-4",
                    H2(:class => "tv-title text-4xl sm:text-5xl", "Four Stages"),
                    P(
                        :class => "tv-body text-lg",
                        "Each stage is one function and one spec struct. A new organ is four more instances, " *
                            "not a new code path."
                    )
                ),
                Div(
                    :class => "grid grid-cols-1 md:grid-cols-2 gap-5 sm:gap-6",
                    [
                        A(
                            :href => "$(BASE)/getting-started/#pipeline",
                            :class => "tv-tile tv-enter", :style => "animation-delay:$(80 * (k - 1))ms",
                            Div(
                                :class => "tv-tile-media aspect-[16/10] flex items-center justify-center p-2 sm:p-3",
                                Img(
                                    :src => "$(BASE)/assets/$(s.image)", :alt => s.alt,
                                    :class => "w-full h-full object-contain rounded-xl"
                                )
                            ),
                            Div(:class => "tv-tile-media px-5 pb-3", P(:class => "tv-tile-caption", s.caption)),
                            Div(
                                :class => "p-6 sm:p-7 space-y-2",
                                Div(
                                    :class => "flex items-baseline gap-3",
                                    Span(:class => "text-sm font-medium tv-accent", s.index),
                                    H3(:class => "tv-title text-2xl sm:text-[1.75rem]", s.title)
                                ),
                                P(:class => "tv-body text-[15px]", s.text)
                            )
                        ) for (k, s) in enumerate(STAGES)
                    ]...
                ),
                Div(
                    :class => "mt-10 flex justify-center",
                    A(:href => "$(BASE)/getting-started/", :class => "tv-pill", "See the Code")
                )
            )
        )
    end
end
