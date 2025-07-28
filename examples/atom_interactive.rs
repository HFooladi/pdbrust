use pdbrust::records::Atom;

fn main() {
    // Create a few test atoms
    let atom1 = Atom::new(
        1,
        "CA".to_string(),
        None,
        "ALA".to_string(),
        "A".to_string(),
        1,
        0.0,
        0.0,
        0.0,
        1.0,
        20.0,
        "C".to_string(),
        None,
    );

    let atom2 = Atom::new(
        2,
        "CA".to_string(),
        None,
        "ALA".to_string(),
        "A".to_string(),
        2,
        1.0,
        0.0,
        0.0,
        1.0,
        20.0,
        "C".to_string(),
        None,
    );

    let atom3 = Atom::new(
        3,
        "CA".to_string(),
        None,
        "ALA".to_string(),
        "A".to_string(),
        3,
        1.0,
        1.0,
        0.0,
        1.0,
        20.0,
        "C".to_string(),
        None,
    );

    let hydrogen = Atom::new(
        4,
        "H".to_string(),
        None,
        "ALA".to_string(),
        "A".to_string(),
        1,
        0.5,
        0.5,
        0.5,
        1.0,
        20.0,
        "H".to_string(),
        None,
    );

    // Print basic atom information
    println!("Atom 1: {:?}", atom1);
    println!("Coordinates: {:?}", atom1.get_coordinates());
    println!("Residue ID: {:?}", atom1.get_residue_id());
    println!("Position: {}", atom1.get_position_string());
    println!("Is backbone: {}", atom1.is_backbone());
    println!("Is hydrogen: {}", atom1.is_hydrogen());
    println!("Is heavy atom: {}", atom1.is_heavy_atom());

    // Distance calculations
    println!("\nDistance between atom1 and atom2: {:.3}", atom1.calculate_distance_to(&atom2));
    println!("Squared distance: {:.3}", atom1.calculate_distance_squared_to(&atom2));

    // Angle calculations
    let angle = atom1.calculate_angle_between(&atom2, &atom3);
    println!("\nAngle between atom1-atom2-atom3: {:.2}°", angle);

    // Different configurations
    let angle2 = atom3.calculate_angle_between(&atom2, &atom1);
    println!("Angle between atom3-atom2-atom1: {:.2}°", angle2);

    println!("\nHydrogen atom: {}", hydrogen.is_hydrogen());
    println!("Distance to hydrogen: {:.3}", atom1.calculate_distance_to(&hydrogen));
} 