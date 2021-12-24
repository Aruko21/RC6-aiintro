#include <string>
#include <iostream>
#include <ctime>
#include <cmath>
#include "SigmoidalNeuron.hpp"
#include "SigmoidalMLP.hpp"

std::ostream &operator<<(std::ostream &os, const std::vector<double> &input) {
    os.precision(3);
    os << "{";
    for (auto const &i: input) {
        os << std::fixed << i << " ";
    }
    os << "}";
    return os;
}

template<typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2){
    std::vector<T> vr(std::begin(v1), std::end(v1));
    vr.insert(std::end(vr), std::begin(v2), std::end(v2));
    return vr;
}

SigmoidalMLPLayer::SigmoidalMLPLayer(size_t neuronsNumber, size_t prevLayerOutputsNum, ConfigurationSingleton& _configuration)
        : configuration(_configuration) {
    for (size_t i = 0; i < neuronsNumber; ++i) {
        this->neurons.emplace_back(SigmoidalNeuron(prevLayerOutputsNum + 1, this->configuration));
    }
};

SigmoidalMLP::SigmoidalMLP(size_t inputSize, size_t outputSize, std::pair<double, double>& _outputScale, ConfigurationSingleton& _configuration)
        : outputScale(_outputScale), configuration(_configuration) {
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

    std::string plotsFilename = configuration.getVariable("PLOTS_FILENAME");

    this->plotsFile.open(plotsFilename);
    if (!this->plotsFile.is_open()) {
        std::cerr << "Error: can't open file " << plotsFilename << std::endl;
        throw std::runtime_error(std::string("Error: can't open file ") + plotsFilename);
    }
}

void SigmoidalMLP::appendPlotToFile(std::vector<std::vector<double>>& plotData) {
    this->plotsFile << std::endl;
    for (size_t point = 0; point < plotData.size(); ++point) {
        for (size_t i = 0; i < plotData.at(point).size(); ++i) {
            this->plotsFile << plotData.at(point).at(i) << " ";
        }
        this->plotsFile << std::endl;
    }
}

void SigmoidalMLP::writePlotDataToFile(const std::string& filename, std::vector<std::vector<double>>& plotData) {
    std::ofstream file(filename);
    for (auto & row : plotData) {
        for (auto & el : row) {
            file << el << " ";
        }
        file << std::endl;
    }
}


double SigmoidalMLP::targetFunction(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& Y) {
    double E = 0.0;
    for (size_t k = 0; k < D.size(); ++k) {
        for (size_t j = 0; j < D.at(k).size(); ++j) {
            E += std::pow(Y.at(k).at(j) - D.at(k).at(j), 2.0);
        }
    }

    return 0.5 * E;
}

double SigmoidalMLP::targetFunctionOnline(std::vector<double>& D, std::vector<double>& Y) {
    double E = 0.0;

    for (size_t j = 0; j < D.size(); ++j) {
        E += std::pow(Y.at(j) - D.at(j), 2.0);
    }

    return 0.5 * E;
}

double SigmoidalMLP::getNormalizedOutput(double output) {
    return (output - this->outputScale.first) * (1 - 0) / (this->outputScale.second - this->outputScale.first);
}
double SigmoidalMLP::getScaledOutput(double output) {
    return output * (this->outputScale.second - this->outputScale.first) / (1 - 0) + this->outputScale.first;
}

void SigmoidalMLP::learn(std::vector<std::vector<double>>& XLearn, std::vector<std::vector<double>>& DLearn) {
    double learningCoef = this->configuration.getVariableDouble("LEARNING_COEF_INIT");
    size_t successIterations = 0;
    double currentTargetFunc = 1.0;

    double learningCoefEPS = this->configuration.getVariableDouble("LEARNING_COEF_EPS");
    double targetFuncEPS = this->configuration.getVariableDouble("TARGET_FUNC_EPS");

//    std::vector<std::vector<double>> normalizedDLearn(DLearn.size());
//    for (size_t k = 0; k < DLearn.size(); ++k) {
//        normalizedDLearn.at(k).resize(DLearn.at(k).size());
//        for (size_t i = 0; i < DLearn.at(k).size(); ++i) {
//            normalizedDLearn.at(k).at(i) = this->getNormalizedOutput(DLearn.at(k).at(i));
//        }
//    }

    // Заполнить файл обучения первоначальным графиком
    std::vector<std::vector<double>> initialMlpPlot(XLearn.size());
    for (size_t k = 0; k < XLearn.size(); ++k) {
        initialMlpPlot.at(k) = XLearn.at(k) + this->getOutput(XLearn.at(k));
    }
    this->appendPlotToFile(initialMlpPlot);
    
    size_t iterations = 1;
    while (learningCoef > learningCoefEPS && currentTargetFunc > targetFuncEPS) {
        // Страшный зверь - выходные значения с прошлого слоя для каждой выборки каждого слоя каждого нейрона
        std::vector<std::vector<std::vector<double>>> prevLayersOutputs;

        // Накопленная ошибка нейрона для каждой выборки каждого нейрона
        std::vector<std::vector<double>> currLayerErrors;
        std::vector<std::vector<double>> leftLayerErrors;

        // Инициализация ошибок на последнем (выходном) слое для каждой выборки
        for (size_t k = 0; k < XLearn.size(); ++k) {
            // Получение результата работы сети для k-й обучающей выборки
            std::vector<double> mlpOutput = this->getOutput(XLearn.at(k));
            prevLayersOutputs.push_back(this->tmpOutputs);
            std::vector<double> errors;

            // Заполнение ошибок (y_i - d_i) для каждого нейрона выходного слоя
            for (size_t i = 0; i < mlpOutput.size(); ++i) {
                errors.push_back(mlpOutput.at(i) - DLearn.at(k).at(i));
            }
            leftLayerErrors.push_back(errors);
        }

        // Массив скорректированных весов для каждого слоя для каждого нейрона
        std::vector<std::vector<std::vector<double>>> correctedWeights(this->layers.size());

        // Для каждого слоя
        for (int l = this->layers.size() - 1; l >= 0; --l) {
            std::vector<SigmoidalNeuron> layerNeurons = this->layers.at(l).getNeurons();
            currLayerErrors = leftLayerErrors;
//            std::cout << "ERRORS for " << l << "layer:" << std::endl;
//            for (auto& errorsForK : currLayerErrors) {
//                std::cout << errorsForK << std::endl;
//            }
//            std::cout << std::endl;

            // Для каждого нейрона в слое
            for (size_t n = 0; n < layerNeurons.size(); ++n) {
                SigmoidalNeuron& neuron = layerNeurons.at(n);
                std::vector<double> weights = neuron.getWeights();
                correctedWeights.at(l).resize(layerNeurons.size());

                if (l != this->layers.size() - 1) {
                    // Умножение ошибок на проивзодную
                    for (size_t k = 0; k < XLearn.size(); ++k) {
                        currLayerErrors.at(k).at(n) *= neuron.getDerivative(prevLayersOutputs.at(k).at(l));
                    }
                }

                // Для каждого синапса (связи)
                for (size_t j = 0; j < neuron.getInputsNumber(); ++j) {
                    correctedWeights.at(l).at(n).resize(neuron.getInputsNumber());

                    double derivative = 0.0;

                    // Для каждой обучающей выборки
                    for (size_t k = 0; k < XLearn.size(); ++k) {
                        // xj - входное значение нейрона связи j
                        double xj = 1.0;
                        // Первая связь - поляризация
                        if (j > 0) {
                            // Выходные значения нейронов на предыдущем слое являются входными для текущего
                            // (j - 1) т.к. j = 0 - поляризация и связи с нейронами начинаются с j = 1
                            xj = prevLayersOutputs.at(k).at(l).at(j - 1);
                        }

                        derivative += currLayerErrors.at(k).at(n) * xj;
                    }
                    // Коррекция веса
                    correctedWeights.at(l).at(n).at(j) = weights.at(j) - learningCoef * derivative;
                }

//                neuron.setWeights(weights);
            }

            // Расчет ошибок для обратного распространения
            if (l > 0) {
                std::vector<SigmoidalNeuron> leftLayerNeurons = this->layers.at(l - 1).getNeurons();

                for (size_t k = 0; k < XLearn.size(); ++k) {
                    leftLayerErrors.at(k).resize(leftLayerNeurons.size());

                    // Для каждого нейрона из предыдущего слоя
                    for (size_t nLeft = 0; nLeft < leftLayerNeurons.size(); ++nLeft) {
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

        std::vector<SigmoidalMLPLayer> correctedLayers;
        for (size_t l = 0; l < this->layers.size(); ++l) {
            size_t numberOfNeurons = this->layers.at(l).getNeurons().size();
            size_t numberOfInputs = this->layers.at(l).getNeurons().at(0).getInputsNumber() - 1;
            correctedLayers.emplace_back(SigmoidalMLPLayer(numberOfNeurons, numberOfInputs, this->configuration));
        }

//        std::copy(this->layers.begin(), this->layers.end(), std::back_inserter(correctedLayers));
        for (size_t l = 0; l < correctedLayers.size(); ++l) {
            std::vector<SigmoidalNeuron>& neurons = correctedLayers.at(l).getNeurons();
            for (size_t n = 0; n < neurons.size(); ++n) {
                SigmoidalNeuron& neuron = neurons.at(n);
                neuron.setWeights(correctedWeights.at(l).at(n));
                SigmoidalNeuron& neuronCompare = neurons.at(n);
            }
        }

        for (size_t k = 0; k < XLearn.size(); ++k) {
            prevE.at(k) = this->getOutput(XLearn.at(k));
            currE.at(k) = this->getOutput(XLearn.at(k), correctedLayers);
        }
        
//        std::vector<std::vector<double>> prevEScaled(DLearn.size());
//        std::vector<std::vector<double>> currEScaled(DLearn.size());

        
//        for (size_t k = 0; k < prevE.size(); ++k) {
//        	prevEScaled.at(k).resize(prevE.at(k).size());
//        	currEScaled.at(k).resize(currE.at(k).size());
//        	for (size_t j = 0; j < prevE.at(k).size(); ++j) {
//        		prevEScaled.at(k).at(j) = this->getScaledOutput(prevE.at(k).at(j));
//        		currEScaled.at(k).at(j) = this->getScaledOutput(currE.at(k).at(j));
//        	}
//        }
                                           
        double targetFuncCurr = this->targetFunction(DLearn, currE);
        double targetFuncPrev = this->targetFunction(DLearn, prevE);

        if (targetFuncCurr < targetFuncPrev) {
            //std::cout << "Successed iteration" << std::endl;
            if (successIterations > 2) {
                learningCoef *= 2.0;
                successIterations = 0;
            } else {
                ++successIterations;
            }

            this->layers = correctedLayers;
            currentTargetFunc = targetFuncCurr;

            /*
            std::cout << "Current learning coefficient: " << learningCoef << std::endl;
            */

            if (iterations % 100 == 0) {
                std::cout << "Current target function value: " << currentTargetFunc << std::endl;
            	std::cout << "EPOCH: " << iterations << std::endl;
            }
            ++iterations;

            std::vector<std::vector<double>> plotData(currE.size());
            for (size_t k = 0; k < currE.size(); ++k) {
                plotData.at(k) = XLearn.at(k) + currE.at(k);
            }
            this->appendPlotToFile(plotData);
        } else {
//            std::cout << "Failed iteration" << std::endl;
            learningCoef /= 2.0;
            successIterations = 0;
        }
    }
    
    std::cout << "Learning completed in " << iterations << " epochs." << std::endl;
}

void SigmoidalMLP::learnOnline(std::vector<double>& XLearn, std::vector<double>& DLearn) {
    double learningCoef = this->configuration.getVariableDouble("LEARNING_COEF_INIT");
    size_t successIterations = 0;
    double currentTargetFunc = 1.0;

    double learningCoefEPS = this->configuration.getVariableDouble("LEARNING_COEF_EPS");
    double targetFuncEPS = this->configuration.getVariableDouble("TARGET_FUNC_EPS");

    while (learningCoef > learningCoefEPS && currentTargetFunc > targetFuncEPS) {
        // выходные значения с прошлого слоя для каждого слоя каждого нейрона
        std::vector<std::vector<double>> prevLayersOutputs;

        // Накопленная ошибка нейрона для каждого нейрона
        std::vector<double> currLayerErrors;
        std::vector<double> leftLayerErrors;

        // Инициализация ошибок на последнем (выходном) слое
        std::vector<double> mlpOutput = this->getOutput(XLearn);
        prevLayersOutputs = this->tmpOutputs;
        std::vector<double> errors;

        // Заполнение ошибок (y_i - d_i) для каждого нейрона выходного слоя
        for (size_t i = 0; i < mlpOutput.size(); ++i) {
            errors.push_back(mlpOutput.at(i) - DLearn.at(i));
        }
        leftLayerErrors = errors;

        // Массив скорректированных весов для каждого слоя для каждого нейрона
        std::vector<std::vector<std::vector<double>>> correctedWeights(this->layers.size());

        // Для каждого слоя
        for (int l = this->layers.size() - 1; l >= 0; --l) {
            std::vector<SigmoidalNeuron> layerNeurons = this->layers.at(l).getNeurons();
            currLayerErrors = leftLayerErrors;
            std::cout << "ERRORS for " << l << "layer: " << currLayerErrors << std::endl;
            std::cout << std::endl;

            // Для каждого нейрона в слое
            for (size_t n = 0; n < layerNeurons.size(); ++n) {
                SigmoidalNeuron& neuron = layerNeurons.at(n);
                std::vector<double> weights = neuron.getWeights();
                correctedWeights.at(l).resize(layerNeurons.size());

                // Умножение ошибок на проивзодную
                for (size_t k = 0; k < XLearn.size(); ++k) {
                    std::cout << "test neuron output: " << neuron.getOutput(prevLayersOutputs.at(l)) << std::endl;
                    std::cout << "test neuron derivative: " <<  neuron.getDerivative(prevLayersOutputs.at(l)) << std::endl;
                    currLayerErrors.at(n) *= neuron.getDerivative(prevLayersOutputs.at(l));
                }

                // Для каждого синапса (связи)
                for (size_t j = 0; j < neuron.getInputsNumber(); ++j) {
                    correctedWeights.at(l).at(n).resize(neuron.getInputsNumber());

//                    double derivative = 0.0;

                    // xj - входное значение нейрона связи j
                    double xj = 1.0;
                    // Первая связь - поляризация
                    if (j > 0) {
                        // Выходные значения нейронов на предыдущем слое являются входными для текущего
                        // (j - 1) т.к. j = 0 - поляризация и связи с нейронами начинаются с j = 1
                        xj = prevLayersOutputs.at(l).at(j - 1);
                    }

                    double derivative = currLayerErrors.at(n) * xj;

                    // Коррекция веса
                    correctedWeights.at(l).at(n).at(j) = weights.at(j) - learningCoef * derivative;
                }
            }

            // Расчет ошибок для обратного распространения
            if (l > 0) {
                std::vector<SigmoidalNeuron> leftLayerNeurons = this->layers.at(l - 1).getNeurons();

                leftLayerErrors.resize(leftLayerNeurons.size());

                // Для каждого нейрона из предыдущего слоя
                for (size_t nLeft = 0; nLeft < leftLayerNeurons.size(); ++nLeft) {
                    double error = 0.0;

                    for (size_t n = 0; n < layerNeurons.size(); ++n) {
                        SigmoidalNeuron& neuron = layerNeurons.at(n);
                        // Взвешиваем ошибку между левым и текущим нейроном
                        error += currLayerErrors.at(n) * neuron.getWeights().at(nLeft + 1);
                    }
                    leftLayerErrors.at(nLeft) = error;
                }

            }
        }

        std::vector<double> prevE(DLearn.size());
        std::vector<double> currE(DLearn.size());

        std::vector<SigmoidalMLPLayer> correctedLayers;
        for (size_t l = 0; l < this->layers.size(); ++l) {
            size_t numberOfNeurons = this->layers.at(l).getNeurons().size();
            size_t numberOfInputs = this->layers.at(l).getNeurons().at(0).getInputsNumber() - 1;
            correctedLayers.emplace_back(SigmoidalMLPLayer(numberOfNeurons, numberOfInputs, this->configuration));
        }

        for (size_t l = 0; l < correctedLayers.size(); ++l) {
            std::vector<SigmoidalNeuron>& neurons = correctedLayers.at(l).getNeurons();
            for (size_t n = 0; n < neurons.size(); ++n) {
                SigmoidalNeuron& neuron = neurons.at(n);
                neuron.setWeights(correctedWeights.at(l).at(n));
                SigmoidalNeuron& neuronCompare = neurons.at(n);
            }
        }


        prevE = this->getOutput(XLearn);
        currE = this->getOutput(XLearn, correctedLayers);

        double targetFuncCurr = this->targetFunctionOnline(DLearn, currE);
        double targetFuncPrev = this->targetFunctionOnline(DLearn, prevE);

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
    std::vector<double> currLayerOutputs;

    this->tmpOutputs.resize(this->layers.size());

    for (size_t l = 0; l < this->layers.size(); ++l) {
        std::vector<SigmoidalNeuron>& layerNeurons = this->layers.at(l).getNeurons();
        currLayerOutputs.resize(layerNeurons.size());
        this->tmpOutputs.at(l).resize(layerNeurons.size());

        if (l != this->layers.size() - 1) {
            for (size_t n = 0; n < layerNeurons.size(); ++n) {
                currLayerOutputs.at(n) = layerNeurons.at(n).getOutput(prevLayerOutputs);
            }
        } else {
            for (size_t n = 0; n < layerNeurons.size(); ++n) {
                auto tmpWeights = layerNeurons.at(n).getWeights();
                currLayerOutputs.at(n) = layerNeurons.at(n).getUSum(prevLayerOutputs, tmpWeights);
            }
        }

        this->tmpOutputs.at(l) = prevLayerOutputs;
        prevLayerOutputs = currLayerOutputs;
    }

    return prevLayerOutputs;
}

std::vector<double> SigmoidalMLP::getOutput(std::vector<double>& X, std::vector<SigmoidalMLPLayer> extLayers) {
    std::vector<double> prevLayerOutputs = X;
    std::vector<double> currLayerOutputs;

    for (size_t l = 0; l < extLayers.size(); ++l) {
        std::vector<SigmoidalNeuron> layerNeurons = extLayers.at(l).getNeurons();
        currLayerOutputs.resize(layerNeurons.size());

        if (l != this->layers.size() - 1) {
            for (size_t n = 0; n < layerNeurons.size(); ++n) {
                currLayerOutputs.at(n) = layerNeurons.at(n).getOutput(prevLayerOutputs);
            }
        } else {
            for (size_t n = 0; n < layerNeurons.size(); ++n) {
                auto tmpWeights = layerNeurons.at(n).getWeights();
                currLayerOutputs.at(n) = layerNeurons.at(n).getUSum(prevLayerOutputs, tmpWeights);
            }
        }

        prevLayerOutputs = currLayerOutputs;
//        currLayerOutputs.clear();
    }

    return prevLayerOutputs;
}

bool SigmoidalMLP::test(std::vector<std::vector<double>>& XTest, std::vector<std::vector<double>>& DTest) {
    std::cout << "testing..." << std::endl;
    double testDiffEPS = this->configuration.getVariableDouble("TEST_DIFF_EPS");

    size_t matchesCount = 0;

    for (size_t i = 0; i < DTest.size(); ++i) {
        std::vector<double> result = this->getOutput(XTest.at(i));
        std::vector<double> scaledResult(result.size());
        for (size_t j = 0; j < result.size(); ++j) {
            scaledResult.at(j) = this->getScaledOutput(result.at(j));
        }
        std::cout << XTest.at(i) << " -> " << scaledResult << " (Real: " << DTest.at(i) << ") |";

        bool matched = true;
        for (size_t j = 0; j < DTest.at(i).size(); ++j) {
            if (std::fabs(scaledResult.at(j) - DTest.at(i).at(j)) >= testDiffEPS) {
                matched = false;
                break;
            }
        }
        if (matched) {
            std::cout << "true";
            ++matchesCount;
        } else {
            std::cout << "false";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Precision of approximation: " << (double) matchesCount / XTest.size() * 100. << "%" << std::endl;

    if (matchesCount != XTest.size()) {
        return false;
    } else {
        return true;
    }
}
