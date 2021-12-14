#ifndef PERCEPTRON_H
#define PERCEPTRON_H

#include <vector>
#include <fstream>

#include "config.hpp"


class SigmoidalNeuron {
private:
    // Количество входов нейрона
    size_t inputsNumber;
    double sigma;
    std::vector<double> weights;
    std::ofstream weightsFile;
    ConfigurationSingleton& configuration;

    double neuronOutput(std::vector<double>& X, std::vector<double>& extWeights);
    double getUSum(std::vector<double>& X, std::vector<double>& extWeights);
    double getDerivative(std::vector<double>& X);
    double targetFunction(std::vector<double>& D, std::vector<double>& Y);
    void appendWeightsToFile();
public:
    explicit SigmoidalNeuron(int _inputs_num, ConfigurationSingleton& _configuration);

    void learn(
        std::vector<std::vector<double>>& XLearn,
        std::vector<double>& DLearn
    );
    bool test(std::vector<std::vector<double>>& XTest, std::vector<double>& DTest);

    std::vector<double> getWeights() { return this->weights; };
    ~SigmoidalNeuron() = default;
};

#endif // PERCEPTRON_H
