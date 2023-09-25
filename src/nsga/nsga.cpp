#include "nsga.hpp"
#include <algorithm>
#include <iostream>
#include <random>
using namespace nsga2;
void NSGA2Solver::initializePopulation() {
  std::uniform_real_distribution<> dis(0.0, 1.0);
  for (int i = 0; i < population.size(); ++i) {
    std::vector<double> individual(dimension);
    std::for_each(individual.begin(), individual.end(),
                  [&, this](auto &x) { x = dis(gen); });
    population[i] = individual;
    auto objectives = objectiveFunction(individual);
    // 评估目标函数值并存储
    if (checkConstraints(individual)) {

      objectiveValues.push_back(objectives);
    } else {
      std::fill(objectives.begin(), objectives.end(), 9999.0);
      // 如果不满足约束，可以考虑如何处理
      objectiveValues.emplace_back(objectives);
    }
  }
}

void NSGA2Solver::evaluatePopulation() {
  // 清空旧的目标函数值

  for (int i = 0; i < population.size(); i++) {
    auto &individual = population[i];
    // 检查约束条件
    if (checkConstraints(individual)) {
      // 计算目标函数值
      objectiveValues[i] = objectiveFunction(individual);

    } else {
      // 如果不满足约束条件，可以考虑如何处理，比如给一个非常差的目标函数值
      objectiveValues.push_back(
          std::vector<double>(objectiveValues[0].size(), 9999.0));
    }
  }
}

void NSGA2Solver::nonDominatedSorting() {
  auto populationSize = population.size();
  fronts.clear();
  std::vector<std::vector<int>> DSet(
      populationSize); // 存储每个个体支配的其他个体
  std::vector<int> n(populationSize, 0); // 存储每个个体被支配的次数
  auto &rank = frontLevels; // 存储每个个体的层级（rank）
  std::fill(rank.begin(), rank.end(), 0);
  // 第一前沿
  std::vector<int> firstFront;

  for (int p = 0; p < populationSize; ++p) {
    const auto &obj_p = objectiveValues[p];
    for (int q = 0; q < populationSize; ++q) {
      if (p == q)
        continue;
      const auto &obj_q = objectiveValues[q];

      // 检查 p 是否支配 q 或 q 是否支配 p
      bool pDominatesQ = true, qDominatesP = true;
      for (size_t i = 0; i < obj_p.size(); ++i) {
        if (obj_p[i] > obj_q[i])
          pDominatesQ = false;
        if (obj_q[i] > obj_p[i])
          qDominatesP = false;
      }

      if (pDominatesQ) {
        DSet[p].push_back(q);
      } else if (qDominatesP) {
        n[p]++;
      }
    }

    if (n[p] == 0) {
      rank[p] = 0;
      firstFront.push_back(p);
    }
  }

  // 初始化前沿
  fronts.push_back(firstFront);

  // 创建其他前沿
  int i = 0;
  while (!fronts[i].empty()) {
    std::vector<int> nextFront;
    for (const int p : fronts[i]) {
      for (const int q : DSet[p]) {
        n[q]--;
        if (n[q] == 0) {
          rank[q] = i + 1;
          nextFront.push_back(q);
        }
      }
    }
    i++;
    fronts.push_back(nextFront);
  }

  // 移除最后一个空的前沿
  if (fronts.back().empty()) {
    fronts.pop_back();
  }
}

int NSGA2Solver::tournamentSelection(int tournamentSize) {
  std::vector<int> selectedCandidates(tournamentSize);
  auto choices = fronts.empty() ? populationSize : fronts[0].size();

  std::uniform_int_distribution<> dist(0, choices - 1); // 均匀分布

  // 从排序后的种群中随机选择 tournamentSize 个候选者
  for (int i = 0; i < tournamentSize; ++i) {
    selectedCandidates[i] =
        populationIndices[dist(gen)]; // 使用类成员 gen 作为随机数生成器
  }

  // 找到这些候选者中最好的一个（即前沿层级最低或拥挤度最高的一个）
  int bestCandidate = selectedCandidates[0];
  for (int i = 1; i < tournamentSize; ++i) {
    int idx1 = selectedCandidates[i];
    int idx2 = bestCandidate;
    if (frontLevels[idx1] < frontLevels[idx2] ||
        (frontLevels[idx1] == frontLevels[idx2] &&
         crowdingDistance[idx1] > crowdingDistance[idx2])) {
      bestCandidate = selectedCandidates[i];
    }
  }

  return bestCandidate;
}
void NSGA2Solver::calculateCrowdingDistance(std::vector<int> &front) {
  // 初始化拥挤度为0
  for (int i : front) {
    crowdingDistance[i] = 0.0;
  }

  int numObjectives = objectiveValues[0].size();

  // 对于每个目标函数进行操作
  for (int m = 0; m < numObjectives; ++m) {
    // 按照第 m 个目标函数对前沿进行排序
    std::sort(front.begin(), front.end(), [&](int i, int j) {
      return objectiveValues[i][m] < objectiveValues[j][m];
    });

    // 设置最端点的拥挤度为无穷大
    crowdingDistance[front[0]] = std::numeric_limits<double>::infinity();
    crowdingDistance[front.back()] = std::numeric_limits<double>::infinity();

    for (size_t i = 1; i < front.size() - 1; ++i) {
      crowdingDistance[front[i]] +=
          (objectiveValues[front[i + 1]][m] - objectiveValues[front[i - 1]][m]);
    }
  }
}
void NSGA2Solver::crossover() {
  // 清空 nextGeneration 以准备新的后代
  nextGeneration.clear();

  // 初始化两个后代
  std::vector<double> offspring1(dimension);
  std::vector<double> offspring2(dimension);
  while (nextGeneration.size() < populationSize) {
    // 通过锦标赛选择两个亲本
    int parent1Index = tournamentSelection(tournamentSize);
    int parent2Index = tournamentSelection(tournamentSize);

    // 获取实际的亲本个体
    const auto &parent1 = population[parent1Index];
    const auto &parent2 = population[parent2Index];

    // 使用 SBX 交叉算子生成后代
    sbxCrossover(parent1, parent2, offspring1, offspring2, eta);

    // 将生成的后代添加到 nextGeneration
    nextGeneration.push_back(offspring1);
    nextGeneration.push_back(offspring2);
  }

  // 如果 nextGeneration 大小超过了预定的种群大小，则进行截断
  if (nextGeneration.size() > populationSize) {
    nextGeneration.resize(populationSize);
  }

  // 在适当的时机，你可以用 nextGeneration 来替换或合并到当前种群 population
}
void NSGA2Solver::mutate() {
  for (auto &individual : nextGeneration) {
    polynomialMutation(individual, m_rate, eta_m); // 例子参数
  }
}
void NSGA2Solver::polynomialMutation(std::vector<double> &individual,
                                     double rate, double eta_m) {
  std::uniform_real_distribution<> dist(0.0, 1.0);

  for (size_t i = 0; i < individual.size(); ++i) {
    if (dist(gen) < rate) { // 以 rate 的概率进行变异
      double u = dist(gen);
      double delta;

      if (u <= 0.5) {
        delta = std::pow(2.0 * u, 1.0 / (eta_m + 1.0)) - 1.0;
      } else {
        delta = 1.0 - std::pow(2.0 * (1.0 - u), 1.0 / (eta_m + 1.0));
      }

      // 应用变异
      individual[i] += delta;
      individual[i] = std::clamp(individual[i], 0.0, 1.0);
      // 注意：这里可能需要根据问题的约束来修正个体
      // 例如，如果某个变量有界限，你可能需要将它裁剪到合法范围内
    }
  }
}
void NSGA2Solver::sbxCrossover(const std::vector<double> &parent1,
                               const std::vector<double> &parent2,
                               std::vector<double> &offspring1,
                               std::vector<double> &offspring2, double eta) {
  std::uniform_real_distribution<> dist(0.0, 1.0);

  for (size_t i = 0; i < parent1.size(); ++i) {
    double u = dist(gen);

    double beta;
    if (u <= 0.5) {
      beta = std::pow(2.0 * u, 1.0 / (eta + 1.0));
    } else {
      beta = std::pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (eta + 1.0));
    }

    // 创建两个后代
    offspring1[i] = 0.5 * ((1 + beta) * parent1[i] + (1 - beta) * parent2[i]);
    offspring2[i] = 0.5 * ((1 - beta) * parent1[i] + (1 + beta) * parent2[i]);
    offspring1[i] = std::clamp(offspring1[i], 0.0, 1.0);
    offspring2[i] = std::clamp(offspring2[i], 0.0, 1.0);
  }
}
void NSGA2Solver::step() {

  crossover();
  mutate();
  evaluatePopulation();
  int nextGenIndex = 0;
  for (int i = populationSize; i < population.size(); i++) {

    population[populationIndices[i]] = nextGeneration[nextGenIndex++];
  }

  nonDominatedSorting();

  for (auto &i : fronts) {
    calculateCrowdingDistance(i);
  }
  int xidx = 0;
  for (int i = 0; i < fronts.size() - 1; ++i) {
    for (int j = 0; j < fronts[i].size(); ++j) {
      populationIndices[xidx++] = fronts[i][j];
    }
  }

  // for (int i = 0; i < populationSize; i++) {
  // auto j = populationIndices[0];
  // std::cout << objectiveValues[j][0] << "," << objectiveValues[j][1]
  //           << std::endl;
  // }
  //对种群进行排序
  std::sort(populationIndices.begin(), populationIndices.end(),
            [&](int i, int j) {
              if (frontLevels[i] == frontLevels[j]) {
                return crowdingDistance[i] > crowdingDistance[j];
              }
              return frontLevels[i] < frontLevels[j];
            });
}
void NSGA2Solver::runNSGA2(int maxGen) {

  initializePopulation();
  for (int i = 0; i < maxGen; i++) {
    step();
  }
}

  std::vector<std::vector<double>> NSGA2Solver::paretoFrontObj() const {

    // 指定切片的起始和终止位置
    auto start = populationIndices.begin(); // 从索引2开始
    auto end = populationIndices.begin() +
               fronts[0].size(); // 到索引5结束（不包括索引6）
    std::vector<std::vector<double>> ret(fronts[0].size());
    int i = 0;
    for (auto it = start; it != end; ++it) {
      ret[i++] = objectiveValues[*it];
    }
    return ret;
  }
  std::vector<std::vector<double>> NSGA2Solver::paretoFrontSol() const {
    // 指定切片的起始和终止位置
    auto start = populationIndices.begin(); // 从索引2开始
    auto end = populationIndices.begin() +
               fronts[0].size(); // 到索引5结束（不包括索引6）
    std::vector<std::vector<double>> ret(fronts[0].size());
    int i = 0;
    for (auto it = start; it != end; ++it) {
      ret[i++] = objectiveValues[*it];
    }
    return ret;
  }