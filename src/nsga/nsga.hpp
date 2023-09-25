#pragma once

#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

namespace nsga2 {

/**
 * @brief Define the type for objective functions.
 */
using ObjectiveFunction = std::function<std::vector<double>(const std::vector<double> &)>;

/**
 * @brief Define the type for constraint functions.
 */
using Constraint = std::function<bool(const std::vector<double> &)>;

class GAParam {
protected:
  double eta = 20.0;      ///< SBX parameter
  int tournamentSize = 3; ///< Tournament size
  double eta_m = 20.0;    ///< SBX parameter
  double m_rate = 0.9;    ///< Mutation rate

public:
  /**
   * @brief Get the value of eta.
   * @return The current value of eta.
   */
  double getEta() const { return eta; }

  /**
   * @brief Set the value of eta.
   * @param newEta The new value of eta to set.
   */
  void setEta(double newEta) { eta = newEta; }

  /**
   * @brief Get the value of tournamentSize.
   * @return The current value of tournamentSize.
   */
  int getTournamentSize() const { return tournamentSize; }

  /**
   * @brief Set the value of tournamentSize.
   * @param newSize The new value of tournamentSize to set.
   */
  void setTournamentSize(int newSize) { tournamentSize = newSize; }

  /**
   * @brief Get the value of eta_m.
   * @return The current value of eta_m.
   */
  double getEtaM() const { return eta_m; }

  /**
   * @brief Set the value of eta_m.
   * @param newEtaM The new value of eta_m to set.
   */
  void setEtaM(double newEtaM) { eta_m = newEtaM; }

  /**
   * @brief Get the value of m_rate.
   * @return The current value of m_rate.
   */
  double getMutationRate() const { return m_rate; }

  /**
   * @brief Set the value of m_rate.
   * @param newMutationRate The new value of m_rate to set.
   */
  void setMutationRate(double newMutationRate) { m_rate = newMutationRate; }
};

class NSGA2Solver : public GAParam {

public:
  /**
   * @brief Perform a single step of the algorithm.
   */
  void step();

  /**
   * @brief Initialize the population.
   */
  void initializePopulation();

  /**
   * @brief Evaluate the population.
   */
  void evaluatePopulation();

  /**
   * @brief Perform non-dominated sorting.
   */
  void nonDominatedSorting();

  /**
   * @brief Calculate crowding distance for individuals.
   * @param indices Indices of the individuals.
   */
  void calculateCrowdingDistance(std::vector<int> &indices);

  /**
   * @brief Select the next generation of individuals.
   */
  void selectNextGeneration();

  /**
   * @brief Perform crossover operation.
   */
  void crossover();

  /**
   * @brief Perform mutation operation.
   */
  void mutate();

  /**
   * @brief Run the NSGA-II algorithm.
   * @param maxGenerations The maximum number of generations to run.
   */
  void runNSGA2(int maxGenerations);

  /**
   * @brief Seed the random number generator.
   * @param seed The seed value for the random number generator.
   */
  void seed(unsigned int seed) {
    gen.seed(seed);
  }

  /**
   * @brief Constructor for NSGA2Solver.
   * @param populationSize The size of the population.
   * @param dimension The dimension of solution vectors.
   */
  NSGA2Solver(int populationSize, int dimension)
      : populationSize(populationSize), dimension(dimension), gen(rd()) {
    // Initialize other members
    crowdingDistance.resize(2 * populationSize, 0.0);
    frontLevels.resize(2 * populationSize, 0);
    population.resize(2 * populationSize);
    populationIndices.resize(2 * populationSize);
    std::iota(populationIndices.begin(), populationIndices.end(),
              0); // Initialize indices
  }

  /**
   * @brief Get the Pareto front objectives.
   * @return A vector of vectors containing Pareto front objectives.
   */
  std::vector<std::vector<double>> paretoFrontObj() const;

  /**
   * @brief Get the Pareto front solutions.
   * @return A vector of vectors containing Pareto front solutions.
   */
  std::vector<std::vector<double>> paretoFrontSol() const;

  ObjectiveFunction objectiveFunction; // Objective function
  Constraint checkConstraints;         // Constraint function
  int populationSize;                  // Population size
  int dimension;                       // Dimension of solution vectors
  std::vector<int> populationIndices;  // Indices of the current population
  std::vector<std::vector<double>> population;      // Current population
  std::vector<std::vector<double>> objectiveValues; // Objective values of the population

private:
  /**
   * @brief Perform tournament selection.
   * @param tournamentSize The size of the tournament.
   * @return The index of the selected individual.
   */
  int tournamentSelection(int tournamentSize);

  /**
   * @brief Perform SBX crossover.
   * @param parent1 The first parent.
   * @param parent2 The second parent.
   * @param offspring1 The first offspring (output).
   * @param offspring2 The second offspring (output).
   * @param eta The SBX parameter.
   */
  void sbxCrossover(const std::vector<double> &parent1,
                    const std::vector<double> &parent2,
                    std::vector<double> &offspring1,
                    std::vector<double> &offspring2, double eta);

  /**
   * @brief Perform polynomial mutation.
   * @param individual The individual to mutate.
   * @param rate The mutation rate.
   * @param eta_m The SBX parameter for mutation.
   */
  void polynomialMutation(std::vector<double> &individual, double rate,
                          double eta_m);

  std::random_device rd;
  std::mt19937 gen;

  std::vector<std::vector<double>> nextGeneration; // Next generation population
  std::vector<double> crowdingDistance;
  std::vector<int> frontLevels;
  std::vector<std::vector<int>> fronts; // Each front is an index list
};

} // namespace nsga2