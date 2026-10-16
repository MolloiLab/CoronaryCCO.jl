# Moved into the package: src/perfusion/flow.jl (Vasculature.Flow). This alias keeps `include(".../flow.jl"); using .Flow`
# working while the scripts migrate; it holds no code.
import Vasculature
const Flow = Vasculature.Flow
