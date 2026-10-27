# NotebookPage — embeds Snapshot's scoped Pluto fragment directly inside the
# Therapy route. Used by per-slug routes registered in docs/app.jl.
#
# `export_notebooks()` writes both a standalone `<slug>.html` and a native
# `<slug>.fragment.html`. The fragment wraps its CSS in `@scope
# (.snap-notebook)`, so Pluto styling cannot leak into the docs shell and the
# notebook reads as ordinary page content instead of a nested iframe viewport.

"""Pretty-print a slug like `00_cardiac_pipeline` → `00 · Cardiac Pipeline`.
Falls back to the raw slug if no override is registered."""
function _notebook_display_title(slug::AbstractString)
    overrides = Dict(
        "00_cardiac_pipeline" => "00 · Cardiac Pipeline",
        "01_pulmonary_pipeline" => "01 · Pulmonary Pipeline",
    )
    return get(overrides, slug, replace(slug, "_" => " "))
end

const _NOTEBOOK_GITHUB_BASE =
    "https://github.com/MolloiLab/Vasculature.jl/blob/main/docs/notebooks"

function NotebookPage(slug::AbstractString)
    BASE = get(ENV, "VASCULATURE_BASE", "")
    title = _notebook_display_title(slug)
    github_src_url = "$(_NOTEBOOK_GITHUB_BASE)/$(slug).jl"
    fragment_path = normpath(joinpath(@__DIR__, "..", "..", "notebooks-static", "$(slug).fragment.html"))
    isfile(fragment_path) || error(
        "Missing notebook fragment $(fragment_path). Run `julia --project=docs/build_env docs/extract_all.jl`."
    )
    notebook_html = replace(
        read(fragment_path, String),
        "__VASCULATURE_NOTEBOOKS_BASE__" => "$(BASE)/notebooks-static",
    )

    return Div(
        :class => "max-w-4xl mx-auto space-y-8",
        # Back pill (left) + source pill (right), then the title.
        Div(
            :class => "flex items-center justify-between gap-3",
            A(:href => "$(BASE)/examples/", :class => "tv-pill tv-pill-sm", "← Examples"),
            A(
                :href => github_src_url, :target => "_blank", :rel => "noopener",
                :class => "tv-pill tv-pill-sm",
                RawHtml(GITHUB_ICON),
                Span(:class => "hidden sm:inline", "View Source"),
            ),
        ),
        H1(:class => "tv-display text-4xl sm:text-5xl", title),

        # Native inline notebook — scoped Pluto/Snapshot styles, one page
        # scroll, inherited Therapy theme, no nested viewport or iframe chrome.
        # The fragment remains real DOM and theme-aware inside a large tile.
        Div(
            :class => "tv-tile px-5 sm:px-10 py-8 sm:py-10",
            RawHtml(notebook_html)
        ),
    )
end
