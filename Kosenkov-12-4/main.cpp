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

template<typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2){
    std::vector<T> vr(std::begin(v1), std::end(v1));
    vr.insert(std::end(vr), std::begin(v2), std::end(v2));
    return vr;
}


int main(int argc, const char** argv) {
    ConfigurationSingleton& configuration = ConfigurationSingleton::getInstance();

    std::vector<std::vector<double>> XLearn;
    std::vector<std::vector<double>> DLearn;

    for (double x = 0; x < 7.0; x += 0.3) {
        XLearn.emplace_back(std::vector<double>{x});
        DLearn.emplace_back(std::vector<double>{std::sin(x)});
    }

    std::string initPlotFilename = configuration.getVariable("INIT_PLOT_FILENAME");
    std::vector<std::vector<double>> initPlotData(XLearn.size());
    for (size_t k = 0; k < XLearn.size(); ++k) {
        initPlotData.at(k) = XLearn.at(k) + DLearn.at(k);
    }
    SigmoidalMLP::writePlotDataToFile(initPlotFilename, initPlotData);

    std::pair<double, double> scale{-1, 1};
    SigmoidalMLP mlpNetwork(1, 1, scale, configuration);

    mlpNetwork.learn(XLearn, DLearn);
    mlpNetwork.test(XLearn, DLearn);


    return 0;
}
