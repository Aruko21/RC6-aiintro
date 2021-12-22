#ifndef SIGMOIDALMLP_H
#define SIGMOIDALMLP_H

#include <vector>
#include <fstream>
#include <cmath>

#include "SigmoidalNeuron.hpp"
#include "config.hpp"

std::ostream &operator<<(std::ostream &os, const std::vector<double> &input);
// Метод для конкатенации двух векторов
template<typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2);


class SigmoidalMLPLayer {
private:
    std::vector<SigmoidalNeuron> neurons;
    ConfigurationSingleton& configuration = ConfigurationSingleton::getInstance();
public:
    SigmoidalMLPLayer() = default;
    SigmoidalMLPLayer& operator=(const SigmoidalMLPLayer& layer) {
        this->neurons = layer.neurons;
        this->configuration = layer.configuration;
        return *this;
    }
    SigmoidalMLPLayer(std::vector<SigmoidalNeuron> _neurons, ConfigurationSingleton& _configuration)
        : neurons(_neurons), configuration(_configuration) {};

    explicit SigmoidalMLPLayer(size_t neuronsNumber, size_t prevLayerNeuronsNum, ConfigurationSingleton& configuration);
    size_t getNumberOfNeurons() { return this->neurons.size(); }
    std::vector<SigmoidalNeuron>& getNeurons() { return this->neurons; };
};

class SigmoidalMLP {
private:
    std::vector<SigmoidalMLPLayer> layers;
    std::ofstream plotsFile;
    std::pair<double, double> outputScale;
    ConfigurationSingleton& configuration;

    double targetFunction(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& Y);
    double targetFunctionOnline(std::vector<double>& D, std::vector<double>& Y);
    double getNormalizedOutput(double output);
    double getScaledOutput(double output);
    void appendPlotToFile(std::vector<std::vector<double>>& plotData);

    // Выходные значения каждого предыдущего слоя каждого нейрона после последнего вызова getOutput
    std::vector<std::vector<double>> tmpOutputs;
public:
    explicit SigmoidalMLP(size_t inputSize, size_t outputSize, std::pair<double, double>& outputScale, ConfigurationSingleton& _configuration);

//    bool test(std::vector<std::vector<double>>& XTest, std::vector<double>& DTest);

    void learn(std::vector<std::vector<double>>& XLearn, std::vector<std::vector<double>>& DLearn);
    void learnOnline(std::vector<double>& XLearn, std::vector<double>& DLearn);
    bool test(std::vector<std::vector<double>>& XTest, std::vector<std::vector<double>>& DTest);
    std::vector<double> getOutput(std::vector<double>& X);
    std::vector<double> getOutput(std::vector<double>& X, std::vector<SigmoidalMLPLayer> extLayers);

    static void writePlotDataToFile(const std::string& filename, std::vector<std::vector<double>>& plotData);

    ~SigmoidalMLP() = default;
};

#endif // SIGMOIDALMLP_H
