#ifndef SIGMOIDALMLP_H
#define SIGMOIDALMLP_H

#include <vector>
#include <fstream>

#include "SigmoidalNeuron.hpp"
#include "config.hpp"


class SigmoidalMLPLayer {
private:
    std::vector<SigmoidalNeuron> neurons;
    ConfigurationSingleton& configuration;
public:
    explicit SigmoidalMLPLayer(size_t neuronsNumber, size_t prevLayerNeuronsNum, ConfigurationSingleton& configuration);
    size_t getNumberOfNeurons() { return this->neurons.size(); }
    std::vector<SigmoidalNeuron> getNeurons() { return this->neurons; };
};

class SigmoidalMLP {
private:
    std::vector<SigmoidalMLPLayer> layers;

    double learnCoef;
    ConfigurationSingleton& configuration;
    double targetFunction(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& Y);

    // Выходные значения каждого предыдущего слоя каждого нейрона после последнего вызова getOutput
    std::vector<std::vector<double>> tmpOutputs;
public:
    explicit SigmoidalMLP(size_t inputSize, size_t outputSize, ConfigurationSingleton& _configuration);

//    bool test(std::vector<std::vector<double>>& XTest, std::vector<double>& DTest);

    void learn(std::vector<std::vector<double>>& XLearn, std::vector<std::vector<double>>& DLearn);
    std::vector<double> SigmoidalMLP::getOutput(std::vector<double>& X);
    std::vector<double> SigmoidalMLP::getOutput(std::vector<double>& X, std::vector<SigmoidalMLPLayer> extLayers);

    ~SigmoidalMLP() = default;
};

#endif // SIGMOIDALMLP_H
