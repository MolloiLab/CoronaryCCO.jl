# Moved into the package: src/perfusion/territory.jl (Vasculature.Territory). This alias keeps `include(".../territory.jl"); using .Territory`
# working while the scripts migrate; it holds no code.
import Vasculature
const Territory = Vasculature.Territory
