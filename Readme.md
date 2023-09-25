# NSGA-II Multi-Objective Optimization

This repository contains a C++ implementation of the NSGA-II algorithm for multi-objective optimization. The algorithm is applied to the ZDT1 problem as a demonstration. Only modern C++ std library is used and no raw pointer is presented.

## Objective Function

The objective function for the ZDT1 problem is implemented in the `calculateZDT1Objectives` function. It takes a vector of decision variables and calculates two objectives.

## Constraint

The `checkConstraints` function is used to check the constraints of the problem. In this example, it ensures that all decision variables are within the range [0, 1].



## Usage

To run the NSGA-II algorithm, follow these steps:

1. Clone this repository.

2. Build the C++ code.

3. Execute the compiled binary.

The Pareto front solutions will be printed to the console.


## License

This code is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- The NSGA-II algorithm is based on the work of K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan in their paper "A Fast Elitist Non-Dominated Sorting Genetic Algorithm for Multi-Objective Optimization: NSGA-II."

Feel free to explore the code and adapt it for your own multi-objective optimization problems.