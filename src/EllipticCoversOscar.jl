module EllipticCoversOscar

using Oscar

export ExplicitMorphism

struct ExplicitMorphism
    source
    target
    x_map
    y_map
    source_coordinates
    target_coordinates
end

formulas(phi::ExplicitMorphism) = (phi.x_map, phi.y_map)

end