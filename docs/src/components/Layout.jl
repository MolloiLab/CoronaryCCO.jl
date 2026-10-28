# Base path read from app.jl (empty in dev, "/Vasculature.jl" in build).
# Hrefs are built as `$(BASE)/path/` so the same Layout works in both modes.
const BASE = get(ENV, "VASCULATURE_BASE", "")

const NAV_ITEMS = [
    ("getting-started", "Getting Started"),
    ("examples", "Examples"),
    ("api", "API"),
]

const GITHUB_ICON = """<svg class="w-[18px] h-[18px]" viewBox="0 0 24 24" fill="currentColor" aria-hidden="true"><path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729 1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997.107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931 0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0 1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.23.653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z"/></svg>"""

"""Vasculature.jl wordmark: the name in the system face, `.jl` dimmed."""
function VasculatureWordmark()
    return NavLink(
        "$(BASE)/",
        RawHtml("""Vasculature<span class="tv-muted">.jl</span>""");
        class = "text-lg font-semibold tracking-tight no-underline hover:opacity-75 transition-opacity",
        active_class = ""
    )
end

function Layout(content)
    return Div(
        :class => "min-h-screen flex flex-col overflow-x-clip transition-colors",
        RawHtml("""<style>html,body{overflow-x:clip}</style>"""),
        # Dark first: with no saved choice the site opens dark (Therapy only
        # adds the class, never removes it). Snapshot's embedded notebook
        # theme follows the same state before the fragment is painted.
        RawHtml("""<script>(function(){var h=document.documentElement;try{var bp=h.getAttribute('data-base-path')||'';if(!localStorage.getItem(bp?'therapy-theme:'+bp:'therapy-theme'))h.classList.add('dark');}catch(e){h.classList.add('dark');}h.setAttribute('data-theme',h.classList.contains('dark')?'fun-dark':'fun-light');})();</script>"""),
        # Sticky top bar: wordmark left, frosted pills right.
        Header(
            :class => "sticky top-0 z-40 h-16 px-4 sm:px-8 backdrop-blur-xl backdrop-saturate-150",
            :style => "background:var(--header-bg)",
            Div(
                :class => "max-w-7xl mx-auto h-full flex items-center justify-between gap-4",
                VasculatureWordmark(),
                Nav(
                    Symbol("aria-label") => "Primary navigation", :class => "flex items-center gap-2",
                    Div(
                        :class => "hidden sm:flex items-center gap-2",
                        [
                            NavLink(
                                "$(BASE)/$(slug)/", label;
                                class = "tv-pill",
                                active_class = "tv-pill-on",
                                inactive_class = ""
                            ) for (slug, label) in NAV_ITEMS
                        ]...
                    ),
                    Details(
                        Symbol("data-mobile-nav") => "", :class => "relative sm:hidden",
                        Summary(
                            Symbol("aria-label") => "Open navigation menu", :class => "tv-pill tv-pill-icon",
                            RawHtml("""<svg class="w-[18px] h-[18px]" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" aria-hidden="true"><path d="M4 7h16M4 12h16M4 17h16"/></svg>""")
                        ),
                        Div(
                            :class => "tv-glass absolute right-0 top-full mt-3 w-52 p-2 flex flex-col gap-1",
                            [
                                NavLink(
                                    "$(BASE)/$(slug)/", label;
                                    class = "tv-pill justify-start w-full",
                                    active_class = "tv-pill-on",
                                    inactive_class = ""
                                ) for (slug, label) in NAV_ITEMS
                            ]...
                        )
                    ),
                    A(
                        :href => "https://github.com/MolloiLab/Vasculature.jl", :target => "_blank",
                        Symbol("aria-label") => "Vasculature.jl on GitHub",
                        :class => "tv-pill tv-pill-icon",
                        RawHtml(GITHUB_ICON)
                    ),
                    DarkModeToggle()
                )
            )
        ),
        # Main content — id="page-content" enables SPA navigation (router swaps this)
        MainEl(
            :id => "page-content", :class => "flex-1 w-full max-w-7xl mx-auto px-4 sm:px-8 pt-10 pb-24 sm:pt-16",
            content
        ),
        Footer(
            :class => "px-4 sm:px-8 py-10",
            Div(
                :class => "max-w-7xl mx-auto flex flex-wrap gap-x-6 gap-y-2 items-center justify-center text-sm tv-muted",
                A(
                    :href => "https://github.com/MolloiLab", :target => "_blank",
                    :class => "no-underline hover:opacity-70 transition-opacity", "MolloiLab"
                ),
                A(
                    :href => "https://github.com/MolloiLab/Vasculature.jl/blob/main/LICENSE", :target => "_blank",
                    :class => "no-underline hover:opacity-70 transition-opacity", "GPL-3.0 License"
                ),
                A(
                    :href => "https://github.com/GroupTherapyOrg/Therapy.jl", :target => "_blank",
                    :class => "no-underline hover:opacity-70 transition-opacity", "Built with Therapy.jl"
                )
            )
        ),
        # Close the compact menu after choosing a destination, clicking away,
        # or pressing Escape. Native <details> keeps keyboard toggling.
        RawHtml("""<script>(function(){var menu=document.querySelector('[data-mobile-nav]');if(!menu)return;document.addEventListener('click',function(e){if(menu.open&&(!menu.contains(e.target)||e.target.closest('[data-mobile-nav] a')))menu.removeAttribute('open');});document.addEventListener('keydown',function(e){if(e.key==='Escape'&&menu.open){menu.removeAttribute('open');var s=menu.querySelector('summary');if(s)s.focus();}});})();</script>"""),
        # Notebook code highlighting, on the first load and after every client-side navigation, of any notebook code
        # block without colour spans. lezer-julia is pinned to the same @lezer/highlight as the highlighter: esm.sh
        # otherwise resolves its ^1.2.0 to a newer copy, whose tags the 1.2.3 classHighlighter does not recognise, and
        # the code stays plain with no error (the page's own highlighter has that fault and marks the code done anyway).
        RawHtml(
            """<script>(function(){var mods=null;function load(){return mods||(mods=Promise.all([import("https://esm.sh/@plutojl/lezer-julia@1.2.0?deps=@lezer/highlight@1.2.3"),import("https://esm.sh/@lezer/highlight@1.2.3")]));}var esc=function(s){return s.replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");};function run(){var todo=Array.prototype.filter.call(document.querySelectorAll("code.pl-jl"),function(el){return !el.querySelector("span");});if(!todo.length)return;load().then(function(m){var parser=m[0].parser,highlightTree=m[1].highlightTree,cls=m[1].classHighlighter;todo.forEach(function(el){if(!el.isConnected||el.querySelector("span"))return;try{var code=el.textContent,tree=parser.parse(code),out="",pos=0;highlightTree(tree,cls,function(from,to,c){if(from>pos)out+=esc(code.slice(pos,from));out+='<span class="'+c+'">'+esc(code.slice(from,to))+"</span>";pos=to;});out+=esc(code.slice(pos));el.innerHTML=out;el.dataset.snapshotHighlighted="1";}catch(e){console.warn("notebook highlight",e);}});}).catch(function(e){console.warn("notebook highlight",e);});}window.addEventListener("therapy:router:loaded",run);document.readyState==="loading"?document.addEventListener("DOMContentLoaded",run):run();})();</script>""",
        ),
    )
end
