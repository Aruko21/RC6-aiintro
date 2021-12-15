/*
 * Вариант 2-7: (Сигмоидальный нейрон)
 *
 * Разработать, программу, моделирующую поведение искусственного трехвходового нейрона указанного типа,
 * и обеспечивающую его обучение для решения задачи классификации.
 *
 * Отладить модель нейрона и процедуру его обучения на произвольных двухмерных данных.
 * Рекомендуется, в тех ситуациях, когда это возможно, использовать режим обучения "оффлайн".
 *
 * Обучить разработанный нейрон на предложенном варианте двухмерных данных и проверить его
 * работу на ряде контрольных точек.
 *
 */

#include <vector>
#include <cmath>
#include <iostream>

#include "config.hpp"
#include "SigmoidalNeuron.hpp"
#include "SigmoidalMLP.hpp"

// Метод для конкатенации двух векторов
template<typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2){
    std::vector<T> vr(std::begin(v1), std::end(v1));
    vr.insert(std::end(vr), std::begin(v2), std::end(v2));
    return vr;
}

//void plot(std::vector<std::vector<double>>& X, SigmoidalNeuron& neuron) {
//    std::ofstream file("plot_f(X).dat");
//    for (auto & row : X) {
//        for (auto & el : row) {
//            file << el << " ";
//        }
//        file << neuron.activation_func(row) << std::endl;
//    }
//}

int main(int argc, const char** argv) {
    ConfigurationSingleton& configuration = ConfigurationSingleton::getInstance();

    std::vector<std::vector<double>> XLearn;
    std::vector<std::vector<double>> DLearn;

    XLearn.emplace_back(std::vector<double>{1.0});
    XLearn.emplace_back(std::vector<double>{2.0});
    XLearn.emplace_back(std::vector<double>{3.0});

    DLearn.emplace_back(std::vector<double>{1.0});
    DLearn.emplace_back(std::vector<double>{2.0});
    DLearn.emplace_back(std::vector<double>{3.0});

//    EllipseDataGenerator::writeToFile(configuration.getVariable("LEARN_DATA_FILENAME"), XLearn, DLearn);
//    EllipseDataGenerator::writeToFile(configuration.getVariable("TEST_DATA_FILENAME"), XTest, DTest);

    SigmoidalMLP mlpNetwork(1, 1, configuration);
//    std::vector<double> test = mlpNetwork.getOutput(XLearn);

//    std::cout << "check output " << test.at(0) << std::endl;
    mlpNetwork.learn(XLearn, DLearn);

//    neuron.test(XTest, DTest);

    return 0;
}
