module EllipticCoversOscar

using Oscar

export ExplicitMorphism, formulas

struct ExplicitMorphism
    source
    target
    x_map
    y_map
    source_coordinates
    target_coordinates
end

function formulas(phi::ExplicitMorphism)
    return phi.x_map, phi.y_map
end

end