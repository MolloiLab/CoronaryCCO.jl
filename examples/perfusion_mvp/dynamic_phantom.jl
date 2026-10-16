# Moved into the package: src/perfusion/dynamic_phantom.jl (Vasculature.DynamicPhantom). This alias keeps `include(".../dynamic_phantom.jl"); using .DynamicPhantom`
# working while the scripts migrate; it holds no code.
import Vasculature
const DynamicPhantom = Vasculature.DynamicPhantom
