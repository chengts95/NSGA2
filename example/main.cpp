#include "nsga.hpp" // Include the necessary header file.
#include <cmath>    // Include the C math library for mathematical functions.
#include <iostream> // Include the standard input/output stream library.
#include <numeric>  // Include the numeric library for accumulate function.
#include <vector> // Include the vector container from the STL (Standard Template Library).

using namespace std; // Use the standard namespace.

// Function to calculate the objectives for the ZDT1 problem.
vector<double> calculateZDT1Objectives(const vector<double> &x) {
  int n = x.size();

  // Objective 1
  double f1 = x[0];

  // Calculate g(x)
  double g = 1.0 + 9.0 * (accumulate(x.begin() + 1, x.end(), 0.0) / (n - 1));

  // Objective 2
  double f2 = g * (1.0 - sqrt(x[0] / g));

  vector<double> objectives = {f1, f2};
  return objectives;
}

// Function to check constraints for a given solution.
bool checkConstraints(const std::vector<double> &solution) {
  for (const auto &x : solution) {
    if (x < 0.0 || x > 1.0) {
      return false;
    }
  }
  return true;
}

// Function to provide a simplified Pareto Front for validation.
vector<vector<double>> paretoFront() {
  // Here, we provide a simplified Pareto Front as an example.
  // In practical applications, you might need more complex Pareto Front data.

  vector<vector<double>> front;

  // Add some solutions to the front.
  front.push_back({0.0, 1.0});
  front.push_back({0.1, 0.9});
  front.push_back({0.2, 0.8});
  front.push_back({0.3, 0.7});
  front.push_back({0.4, 0.6});
  front.push_back({0.5, 0.5});

  return front;
}

int main() {
  // Create an instance of the NSGA2Solver with a population size of 100 and 30
  // dimensions.
  nsga2::NSGA2Solver solver(200, 30);

  // Seed the random number generator with 100.
  solver.seed(1);


  // Set the objective function to calculateZDT1Objectives.
  solver.objectiveFunction = calculateZDT1Objectives;

  // Set the constraint-checking function to checkConstraints.
  solver.checkConstraints = checkConstraints;

  // Run the NSGA-II algorithm.
  solver.runNSGA2(40);

  // Get the Pareto front solutions.
  auto sol = solver.paretoFrontObj();

  // Print the Pareto front solutions.
  for (int i = 0; i < sol.size(); i++) {
    cout << sol[i][0] << "," << sol[i][1] << endl;
  }

  return 0;
}