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
    ConfigurationSingleton& configuration;

    double neuronOutput(std::vector<double>& X, std::vector<double>& extWeights);

    double targetFunction(std::vector<double>& D, std::vector<double>& Y);
public:
    explicit SigmoidalNeuron(int _inputs_num, ConfigurationSingleton& _configuration);
    SigmoidalNeuron& operator=(const SigmoidalNeuron& neuron) {
        this->inputsNumber = neuron.inputsNumber;
        this->sigma = neuron.sigma;
        this->weights = neuron.weights;
        this->configuration = neuron.configuration;
        return *this;
    };


    static double generateRandNumInRange(double l, double r);

    void learn(
        std::vector<std::vector<double>>& XLearn,
        std::vector<double>& DLearn
    );
    double getUSum(std::vector<double>& X, std::vector<double>& extWeights);
    bool test(std::vector<std::vector<double>>& XTest, std::vector<double>& DTest);
    double getOutput(std::vector<double>& X);
    double getDerivative(std::vector<double>& X);
    size_t getInputsNumber() { return this->inputsNumber; };
    std::vector<double> getWeights() { return this->weights; };
    void setWeights(std::vector<double>& newWeights) { this->weights = newWeights; };
    ~SigmoidalNeuron() = default;
};

#endif // PERCEPTRON_H
