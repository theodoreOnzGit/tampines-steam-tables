extern crate ndarray;
extern crate ndarray_linalg;

use ndarray::{array, Array1, Array2};
use ndarray_linalg::Solve; // For potential solving, though not directly applicable here due to singularity

#[test]
fn steam_mass_bal_vibe_code() {
    // Define the coefficient matrix A for mass balance
    // The rows correspond to the control volumes (Boiler, Turbine, Condenser, Pump)
    // The columns correspond to the mass flow rates (m1, m2, m3, m4)
    let a: Array2<f64> = array![
        [1.0, -1.0, 0.0, 0.0], // Boiler: m1 - m2 = 0
        [0.0, 1.0, -1.0, 0.0],  // Turbine: m2 - m3 = 0
        [0.0, 0.0, 1.0, -1.0],  // Condenser: m3 - m4 = 0
        [-1.0, 0.0, 0.0, 1.0],  // Pump: m4 - m1 = 0
    ];

    // Define the right-hand side vector b (all zeros for simple mass balance)
    let b: Array1<f64> = array![0.0, 0.0, 0.0, 0.0];

    println!("Mass Balance Coefficient Matrix (A):\n{}", a);
    println!("\nRight-hand side vector (b):\n{}", b);

    // --- Discussion on Solving ---
    println!("\n--- Attempting to solve A * x = b ---");
    // For a singular matrix like this, directly calling a solve method
    // (e.g., a.solve(&b)) will typically fail or indicate singularity.
    // Let's check the determinant to confirm singularity.
    match a.det() {
        Ok(det_val) => {
            println!("Determinant of A: {}", det_val);
            if det_val.abs() < 1e-9 { // Check if determinant is close to zero
                println!("Matrix A is singular (determinant is approximately zero).");
                println!("This means the system of equations A * x = 0 has infinitely many solutions.");
                println!("Physically, this implies that the mass flow rate is constant throughout the cycle:");
                println!("m1 = m2 = m3 = m4 = k, where k can be any arbitrary mass flow rate.");
                println!("To get a specific numerical value for the mass flow rate, you would need to provide an external constraint,");
                println!("such as specifying one of the mass flow rates (e.g., m1 = 10.0 kg/s).");
                println!("If you were to add such a constraint, you would modify the matrix and vector accordingly.");
            } else {
                println!("Matrix A is non-singular.");
                // If it were non-singular, you could solve it like this:
                // match a.solve(&b) {
                //     Ok(x_solution) => println!("\nSolution vector (x):\n{}", x_solution),
                //     Err(e) => println!("\nError solving system: {:?}", e),
                // }
            }
        },
        Err(e) => println!("Error calculating determinant: {:?}", e),
    }

    // Example of how you might modify the system if you *knew* m1, for instance.
    // If m1 is known, say m1 = 10.0, then the first equation becomes a constraint.
    // You could then reduce the system or modify the matrix.
    // For example, if m1 is known, you could substitute it into the last equation,
    // and then you'd have a 3x3 system for m2, m3, m4.
    println!("\n--- Example with an additional constraint (m1 = 10.0 kg/s) ---");
    // Let's assume m1 = 10.0.
    // We can modify the first equation to be m1 = 10.0
    // And the last equation to be -m1 + m4 = 0 => m4 = m1
    // This is one way to handle it, by making one equation a direct assignment.
    // A more formal way is to reduce the system.
    // For simplicity, let's just state the consequence:
    println!("If we set m1 = 10.0 kg/s, then due to mass balance:");
    println!("m2 = 10.0 kg/s (from boiler)");
    println!("m3 = 10.0 kg/s (from turbine)");
    println!("m4 = 10.0 kg/s (from condenser)");
    println!("And this satisfies the pump equation (m4 - m1 = 0).");
}
