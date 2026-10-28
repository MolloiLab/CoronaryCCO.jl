# TableOfContents.jl — right-hand list of the page's sections.
#
# A column of its own at xl+ breakpoints, sticky under the top bar and aligned with the bar's right edge, so it never
# floats over the content at any width; hidden on smaller viewports. Each entry is an anchor link to a section `id`
# within the current page.

function TableOfContents(sections::Vector{Tuple{String, String}})
    return Nav(
        :class => "hidden xl:block sticky top-28 self-start max-h-[calc(100vh-8rem)] overflow-y-auto",
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
    return Div(
        :class => "xl:grid xl:grid-cols-[minmax(0,1fr)_12rem] xl:gap-16",
        Div(:class => "min-w-0", content),
        TableOfContents(sections),
    )
end
