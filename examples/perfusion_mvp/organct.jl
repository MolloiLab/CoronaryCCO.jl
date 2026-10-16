# Moved into the package: src/perfusion/organct.jl (Vasculature.OrganCT). This alias keeps `include(".../organct.jl"); using .OrganCT`
# working while the scripts migrate; it holds no code.
import Vasculature
const OrganCT = Vasculature.OrganCT
