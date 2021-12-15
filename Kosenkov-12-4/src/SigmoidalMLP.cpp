#include <string>
#include <iostream>
#include <ctime>
#include "SigmoidalNeuron.hpp"
#include "SigmoidalMLP.hpp"


SigmoidalMLPLayer::SigmoidalMLPLayer(size_t neuronsNumber, size_t prevLayerOutputsNum, ConfigurationSingleton& _configuration)
        : configuration(_configuration) {
    for (size_t i = 0; i < neuronsNumber; ++i) {
        this->neurons.emplace_back(SigmoidalNeuron(prevLayerOutputsNum + 1, this->configuration));
    }
};

SigmoidalMLP::SigmoidalMLP(size_t inputSize, size_t outputSize, ConfigurationSingleton& _configuration)
        : configuration(_configuration) {
    srand(time(nullptr));
    size_t layersNum = configuration.getVariableInt("MLP_LAYERS_NUMBER");

    for (size_t i = 0; i < layersNum; ++i) {
        size_t neuronsNum = configuration.getVariableInt("MLP_LAYER" + std::to_string(i + 1) + "_NEURONS");
        // Если слой первый - то количество входов соотвтвует размерности данных
        size_t prevLayerOutputsNum = inputSize;
        if (i > 0) {
            // Иначе - количеству нейронов на прошлом слое
            prevLayerOutputsNum = configuration.getVariableInt("MLP_LAYER" + std::to_string(i) + "_NEURONS");
        }
        this->layers.emplace_back(SigmoidalMLPLayer(neuronsNum, prevLayerOutputsNum, this->configuration));
    }

    size_t lastLayerSize = this->layers.at(this->layers.size() - 1).getNumberOfNeurons();
    this->layers.emplace_back(SigmoidalMLPLayer(outputSize, lastLayerSize, this->configuration));
}

double SigmoidalMLP::targetFunction(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& Y) {
    double E = 0.0;
    for (size_t k = 0; k < D.size(); ++k) {
        for (size_t j = 0; D.at(k).size(); ++j) {
            E += std::pow(Y.at(k).at(j) - D.at(k).at(j), 2.0);
        }
    }

    return 0.5 * E;
}

void SigmoidalMLP::learn(std::vector<std::vector<double>>& XLearn, std::vector<std::vector<double>>& DLearn) {
    double learningCoef = this->configuration.getVariableDouble("LEARNING_COEF_INIT");
    size_t successIterations = 0;
    double currentTargetFunc = 1.0;

    double learningCoefEPS = this->configuration.getVariableDouble("LEARNING_COEF_EPS");
    double targetFuncEPS = this->configuration.getVariableDouble("TARGET_FUNC_EPS");

    while (learningCoef > learningCoefEPS && currentTargetFunc > targetFuncEPS) {
        // Страшный зверь - выходные значения для каждой выборки каждого слоя каждого нейрона
        std::vector<std::vector<std::vector<double>>> prevLayersOutputs;

        // Накопленная ошибка нейрона для каждой выборки каждого нейрона
        std::vector<std::vector<double>> currLayerErrors;
        std::vector<std::vector<double>> leftLayerErrors;

        for (size_t k = 0; k < XLearn.size(); ++k) {
            std::vector<double> mlpOutput = this->getOutput(XLearn.at(k));
            prevLayersOutputs.push_back(this->tmpOutputs);
            std::vector<double> errors;

            // Заполнение ошибок (y_i - d_i) для каждого нейрона выходного слоя
            for (size_t i = 0; i < mlpOutput.size(); ++i) {
                errors.push_back(mlpOutput.at(i) - DLearn.at(k).at(i));
            }
            currLayerErrors.push_back(errors);
        }

        std::vector<std::vector<std::vector<double>>> correctedWeights(this->layers.size());

        // Для каждого слоя
        for (size_t l = this->layers.size() - 1; l >= 0; --l) {
            std::vector<SigmoidalNeuron> layerNeurons = this->layers.at(l).getNeurons();
            currLayerErrors = leftLayerErrors;

            // Для каждого нейрона в слое
            for (size_t n = 0; n < layerNeurons.size(); ++n) {
                SigmoidalNeuron& neuron = layerNeurons.at(n);
                std::vector<double> weights = neuron.getWeights();
                correctedWeights.at(l).resize(layerNeurons.size());

                // Для каждого синапса
                for (size_t j = 0; j < neuron.getInputsNumber(); ++j) {
                    correctedWeights.at(l).at(n).resize(neuron.getInputsNumber());

                    double derivative = 0.0;

                    // Для каждой обучающей выборки
                    for (size_t k = 0; k < XLearn.size(); ++k) {
                        double xj = 1.0;
                        if (j > 0) {
                            // Выходные значения нейронов на предыдущем слое являются входными для текущего
                            xj = prevLayersOutputs.at(k).at(l).at(j - 1);
                        }

                        // Берем дельту для нейрона
                        currLayerErrors.at(k).at(n) *= neuron.getDerivative(XLearn.at(k));
                        derivative += currLayerErrors.at(k).at(n) * xj;
                    }
                    // Коррекция веса
                    correctedWeights.at(l).at(n).at(j) = weights.at(j) - learningCoef * derivative;
                }

                neuron.setWeights(weights);
            }


            if (l > 0) {
                std::vector<SigmoidalNeuron> leftLayerNeurons = this->layers.at(l - 1).getNeurons();
                leftLayerErrors.resize(leftLayerNeurons.size());

                for (size_t nLeft = 0; nLeft < leftLayerNeurons.size(); ++nLeft) {
                    for (size_t k = 0; k < XLearn.size(); ++k) {
                        double error = 0.0;

                        for (size_t n = 0; n < layerNeurons.size(); ++n) {
                            SigmoidalNeuron& neuron = layerNeurons.at(n);
                            // Взвешиваем ошибку синапсом между левым и текущим нейроном
                            error += currLayerErrors.at(k).at(n) * neuron.getWeights().at(nLeft + 1);
                        }
                        leftLayerErrors.at(k).at(nLeft) = error;
                    }
                }
            }
        }

        std::vector<std::vector<double>> prevE(DLearn.size());
        std::vector<std::vector<double>> currE(DLearn.size());

        std::vector<SigmoidalMLPLayer> correctedLayers = this->layers;
        for (size_t l = 0; l < correctedLayers.size(); ++l) {
            std::vector<SigmoidalNeuron> neurons = correctedLayers.at(l).getNeurons();
            for (size_t n = 0; n < neurons.size(); ++n) {
                SigmoidalNeuron neuron = neurons.at(n);
                neuron.setWeights(correctedWeights.at(l).at(n));
            }
        }

        for (size_t k = 0; k < XLearn.size(); ++k) {
            prevE.at(k) = this->getOutput(XLearn.at(k));
            currE.at(k) = this->getOutput(XLearn.at(k), correctedLayers);
        }

        double targetFuncCurr = this->targetFunction(DLearn, currE);
        double targetFuncPrev = this->targetFunction(DLearn, prevE);

        if (targetFuncCurr < targetFuncPrev) {
            std::cout << "Successed iteration" << std::endl;
            if (successIterations > 2) {
                learningCoef *= 2.0;
                successIterations = 0;
            } else {
                ++successIterations;
            }

            this->layers = correctedLayers;
            currentTargetFunc = targetFuncCurr;

            std::cout << "Current learning coefficient: " << learningCoef << std::endl;
            std::cout << "Current target function value: " << currentTargetFunc << std::endl;
        } else {
            std::cout << "Failed iteration" << std::endl;
            learningCoef /= 2.0;
            successIterations = 0;
        }
    }
}

// Возвращает выходные значения для каждого нейрона в выходном слое
std::vector<double> SigmoidalMLP::getOutput(std::vector<double>& X) {
    std::vector<double> prevLayerOutputs = X;
    std::vector<double> nextLayerOutputs;

    this->tmpOutputs.clear();

    for (size_t i = 0; i < this->layers.size(); ++i) {
        std::vector<SigmoidalNeuron>& layerNeurons = this->layers.at(i).getNeurons();
        for (size_t n = 0; n < layerNeurons.size(); ++n) {
            nextLayerOutputs.push_back(layerNeurons.at(n).getOutput(prevLayerOutputs));
        }

        this->tmpOutputs.push_back(prevLayerOutputs);
        prevLayerOutputs = nextLayerOutputs;
        nextLayerOutputs.clear();
    }

    return prevLayerOutputs;
}

std::vector<double> SigmoidalMLP::getOutput(std::vector<double>& X, std::vector<SigmoidalMLPLayer> extLayers) {
    std::vector<double> prevLayerOutputs = X;
    std::vector<double> nextLayerOutputs;

    for (size_t i = 0; i < extLayers.size(); ++i) {
        std::vector<SigmoidalNeuron>& layerNeurons = extLayers.at(i).getNeurons();
        for (size_t n = 0; n < layerNeurons.size(); ++n) {
            nextLayerOutputs.push_back(layerNeurons.at(n).getOutput(prevLayerOutputs));
        }

        prevLayerOutputs = nextLayerOutputs;
        nextLayerOutputs.clear();
    }

    return prevLayerOutputs;
}