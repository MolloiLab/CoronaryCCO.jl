# Moved into the package: src/perfusion/body.jl (Vasculature.Body). This alias keeps `include(".../body.jl"); using .Body`
# working while the scripts migrate; it holds no code.
import Vasculature
const Body = Vasculature.Body
