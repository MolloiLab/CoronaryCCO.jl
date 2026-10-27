# TableOfContents.jl — floating right-hand list of the page's sections.
#
# Fixed to the right side at xl+ breakpoints; hidden on smaller viewports.
# Each entry is an anchor link to a section `id` within the current page.

function TableOfContents(sections::Vector{Tuple{String, String}})
    return Nav(
        :class => "hidden xl:block fixed right-10 top-28 w-48",
        Div(
            :class => "space-y-1",
            P(:class => "text-xs font-medium tv-muted mb-3 px-3", "On This Page"),
            For(sections) do (id, label)
                A(
                    :href => "#$id",
                    :class => "block rounded-full px-3 py-1.5 text-[13px] tv-body no-underline hover:bg-[var(--pill-hover)] transition-colors",
                    label
                )
            end
        )
    )
end

function PageWithTOC(sections::Vector{Tuple{String, String}}, content)
    return Fragment(
        content,
        TableOfContents(sections),
    )
end
