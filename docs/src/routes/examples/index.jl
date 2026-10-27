# Examples gallery — auto-discovers Pluto notebooks from docs/notebooks/.
# Each tile links to /examples/<slug>/, which is registered programmatically
# in docs/app.jl and served via NotebookPage(slug).
#
# Tile metadata (title, summary, thumbnail) lives in `NOTEBOOK_META` below. A
# notebook without an entry still gets a tile built from its slug.

let BASE = get(ENV, "VASCULATURE_BASE", ""),
        notebooks_dir = joinpath(@__DIR__, "..", "..", "..", "notebooks"),
        NOTEBOOK_META = Dict(
        "00_cardiac_pipeline" => (
            index = "00",
            title = "Cardiac Pipeline",
            summary = "From the raw XCAT heart to a perfused CT phantom, and the same heart with a 90% LAD stenosis.",
            thumbnail = "step4_heart_phantom.png",
            tags = ["Heart", "XCAT", "Phantom"],
        ),
    )

    # Tile builder — defined inside the `let` so Therapy's file-based router
    # doesn't try to register it as a route.
    notebook_tile = function (slug::AbstractString, k::Integer)
        meta = get(
            NOTEBOOK_META, slug, (
                index = first(split(slug, "_")),
                title = replace(slug, "_" => " "),
                summary = "",
                thumbnail = nothing,
                tags = String[],
            )
        )

        return A(
            :href => "$(BASE)/examples/$(slug)/",
            :class => "tv-tile tv-enter", :style => "animation-delay:$(80 * (k - 1))ms",
            if meta.thumbnail !== nothing
                Div(
                    :class => "tv-tile-media aspect-[16/9] flex items-center justify-center p-4",
                    Img(
                        :src => "$(BASE)/assets/$(meta.thumbnail)", :alt => meta.title,
                        :class => "max-w-full max-h-full object-contain rounded-xl"
                    )
                )
            else
                Div(:class => "aspect-[16/9]", :style => "background:var(--surface-2)")
            end,
            Div(
                :class => "p-6 sm:p-7 space-y-3",
                Div(
                    :class => "flex items-baseline gap-3",
                    Span(:class => "text-sm font-medium tv-accent", meta.index),
                    H3(:class => "tv-title text-2xl", meta.title)
                ),
                P(:class => "tv-body text-[15px]", meta.summary),
                Div(
                    :class => "flex flex-wrap items-center gap-2 pt-1",
                    [Span(:class => "tv-pill tv-pill-sm pointer-events-none", tag) for tag in meta.tags]...
                ),
                Div(:class => "tv-tile-more pt-2 text-sm font-medium", "Open Notebook →")
            )
        )
    end

    () -> begin
        # Discover slugs at render time so a fresh notebook drop shows up
        # without code changes here. Like app.jl, leave out notebooks that are
        # skipped or not rendered yet, so no tile points at a missing page.
        skip = Set(strip.(split(get(ENV, "VASCULATURE_SKIP_NOTEBOOKS", ""), ","; keepempty = false)))
        rendered(s) = isfile(joinpath(notebooks_dir, "..", "notebooks-static", "$(s).fragment.html"))
        slugs = isdir(notebooks_dir) ?
            sort(
                [
                    splitext(f)[1] for f in readdir(notebooks_dir)
                    if endswith(f, ".jl") && !endswith(f, ".sessions.toml") && !(splitext(f)[1] in skip) &&
                    rendered(splitext(f)[1])
                ]
            ) :
            String[]

        Div(
            :class => "space-y-12",
            Div(
                :class => "max-w-2xl space-y-4",
                H1(:class => "tv-display text-5xl sm:text-6xl", "Examples"),
                P(
                    :class => "tv-body text-lg",
                    "Worked examples of the Vasculature.jl pipeline. Each one is a Pluto notebook, shown with its ",
                    "code, text and figures in place."
                )
            ),
            if isempty(slugs)
                Div(
                    :class => "tv-tile py-20 text-center",
                    P(
                        :class => "tv-body",
                        "No rendered notebooks yet. Add a .jl file to ",
                        Code(:class => "tv-code", "docs/notebooks/"),
                        " and rebuild."
                    )
                )
            else
                Div(
                    :class => "grid grid-cols-1 md:grid-cols-2 gap-5 sm:gap-6",
                    [notebook_tile(slug, k) for (k, slug) in enumerate(slugs)]...
                )
            end
        )
    end
end
