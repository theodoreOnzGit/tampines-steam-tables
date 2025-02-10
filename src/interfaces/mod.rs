/// this set of interfaces allows the user 
/// to interact using a more functional programming 
/// style (no objects)
///
/// this keeps things simple.
pub mod functional_programming;

/// for OOP users who want to make a struct (class)
/// and then use that for extracting data, 
/// this is where the stuff is stored
pub mod object_oriented_programming;

/// these tests show you how to use the interfaces 
///
/// i may attempt to produce part or whole of the 
/// steam tables here
#[cfg(test)]
pub mod tests_and_examples;
