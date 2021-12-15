#include <string>
#include <iostream>
#include <chrono>
#include <cmath>
#include "SigmoidalNeuron.hpp"

double SigmoidalNeuron::generateRandNumInRange(double l, double r) {
    double tmp = (double) std::rand() / (RAND_MAX / (r - l));
    return l + tmp;
}

SigmoidalNeuron::SigmoidalNeuron(int _inputsNumber, ConfigurationSingleton& _configuration)
        : inputsNumber(_inputsNumber), configuration(_configuration) {
    this->sigma = configuration.getVariableDouble("SIGMA");

    std::string weightsFilename = configuration.getVariable("WEIGHTS_FILENAME");

    this->weights.resize(inputsNumber);

    for (size_t i = 0; i < inputsNumber; ++i) {
        weights.at(i) = SigmoidalNeuron::generateRandNumInRange(0, 1);
    }
}

double SigmoidalNeuron::neuronOutput(std::vector<double>& X, std::vector<double>& extWeights) {
    // Функция активации
    return 1 / ( 1 + std::exp(- this->sigma * this->getUSum(X, extWeights)) );
}

double SigmoidalNeuron::getOutput(std::vector<double>& X) {
    return this->neuronOutput(X, this->weights);
}

double SigmoidalNeuron::getUSum(std::vector<double>& X, std::vector<double>& extWeights) {
    double u = 0;
    // X0 - вход поляризации = 1
    u += extWeights.at(0) * 1.0;

    for (size_t i = 1; i < this->inputsNumber; ++i) {
        u += X.at(i - 1) * extWeights.at(i);
    }

    return u;
}

double SigmoidalNeuron::getDerivative(std::vector<double>& X) {
    double output = this->neuronOutput(X, this->weights);
    return this->sigma * output * (1 - output);
}

double SigmoidalNeuron::targetFunction(std::vector<double>& D, std::vector<double>& Y) {
    double E = 0.0;
    for (size_t k = 0; k < D.size(); ++k) {
        E += std::pow(Y.at(k) - D.at(k), 2.0);
    }

    return 0.5 * E;
}

void SigmoidalNeuron::learn(
    std::vector<std::vector<double>>& XLearn,
    std::vector<double>& DLearn
) {
    std::cout << "learning..." << std::endl;
    std::vector<double> correctedWeights(this->weights.size());

    double learningCoef = this->configuration.getVariableDouble("LEARNING_COEF_INIT");
    size_t successIterations = 0;
    double currentTargetFunc = 1.0;

    double learningCoefEPS = this->configuration.getVariableDouble("LEARNING_COEF_EPS");
    double targetFuncEPS = this->configuration.getVariableDouble("TARGET_FUNC_EPS");

    while (learningCoef > learningCoefEPS && currentTargetFunc > targetFuncEPS) {
        for (size_t j = 0; j < this->inputsNumber; ++j) {
            double derivative = 0.0;

            for (size_t k = 0; k < XLearn.size(); ++k) {
                double neuronOut = this->neuronOutput(XLearn.at(k), this->weights);

                double xj = 1.0;
                if (j > 0) {
                    xj = XLearn.at(k).at(j - 1);
                }
                // Вычисление компоненты градиента
                derivative += (neuronOut - DLearn.at(k)) * this->getDerivative(XLearn.at(k)) * xj;
            }

            // коорекция веса
            correctedWeights.at(j) = this->weights.at(j) - learningCoef * derivative;
        }

        std::vector<double> prevE(DLearn.size());
        std::vector<double> currE(DLearn.size());

        for (size_t k = 0; k < XLearn.size(); ++k) {
            prevE.at(k) = this->neuronOutput(XLearn.at(k), this->weights);
            currE.at(k) = this->neuronOutput(XLearn.at(k), correctedWeights);
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

            this->weights = correctedWeights;
            currentTargetFunc = targetFuncCurr;

            std::cout << "Corrected weights: " << std::endl;
            std::cout << "Current learning coefficient: " << learningCoef << std::endl;
            std::cout << "Current target function value: " << currentTargetFunc << std::endl;
        } else {
            std::cout << "Failed iteration" << std::endl;
            learningCoef /= 2.0;
            successIterations = 0;
        }
    }
}

bool SigmoidalNeuron::test(std::vector<std::vector<double>>& XTest, std::vector<double>& DTest) {
    std::cout << "testing..." << std::endl;
    double testDiffEPS = this->configuration.getVariableDouble("TEST_DIFF_EPS");

    size_t matchesCount = 0;

    for (size_t i = 0; i < DTest.size(); ++i) {
        double result = this->neuronOutput(XTest.at(i), this->weights);
        std::cout << result << " = " << DTest.at(i) << "? -> ";

        if (std::fabs(result -  DTest.at(i)) < testDiffEPS) {
            std::cout << "true";
            ++matchesCount;
        } else {
            std::cout << "false";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Precision of classification: " << (double) matchesCount / XTest.size() * 100. << "%" << std::endl;

    if (matchesCount != XTest.size()) {
        return false;
    } else {
        return true;
    }
}
