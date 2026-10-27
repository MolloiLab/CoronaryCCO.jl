# CodeBlock — a dark code listing with a copy badge in its top-right corner (the pattern of BasisSimulator.jl's docs).
# The badge copies the block's code and marks itself `data-copied` for a moment; CSS (`.tv-copy`) swaps its label.

const _COPY_BUTTON = """<button type="button" aria-label="Copy to the clipboard" onclick="(function(b){var c=b.closest('[data-codeblock]').querySelector('code');navigator.clipboard.writeText(c.innerText).then(function(){b.setAttribute('data-copied','');setTimeout(function(){b.removeAttribute('data-copied')},1400)})})(this)" class="tv-pill tv-pill-sm tv-copy"><span class="tv-copy-idle">Copy</span><span class="tv-copy-done">Copied</span></button>"""

"""
    CodeBlock(code; language = "julia", id = nothing)

A dark code listing with a copy badge. `id` names the listing's `<pre>`, for links to it.
"""
function CodeBlock(code::AbstractString; language::AbstractString = "julia", id = nothing)
    return Div(
        Symbol("data-codeblock") => "", :class => "relative",
        RawHtml(_COPY_BUTTON),
        id === nothing ? Pre(:class => "tv-pre", Code(:class => "language-$language", code)) :
            Pre(:class => "tv-pre", :id => id, Code(:class => "language-$language", code)),
    )
end
