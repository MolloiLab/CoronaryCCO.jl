# Moved into the package: src/perfusion/cfc.jl (Vasculature.CFC). This alias keeps `include(".../cfc.jl"); using .CFC`
# working while the scripts migrate; it holds no code.
import Vasculature
const CFC = Vasculature.CFC
