# Moved into the package: Vasculature.couple_veins (src/perfusion/pipeline.jl). This alias keeps
# `include(".../coupling.jl")` working while the scripts migrate; it holds no code.
import Vasculature
const couple_veins = Vasculature.couple_veins
