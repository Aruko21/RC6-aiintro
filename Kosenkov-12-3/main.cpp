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
#include "EllipseDataGenerator.hpp"
#include "EllipseParams.hpp"

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
    EllipseParams ellipse1(configuration, "1");
    EllipseParams ellipse2(configuration, "2");

    // Генерация данных в эллипсах для обучения
    EllipseDataGenerator ellipse1Generator(ellipse1);
    EllipseDataGenerator ellipse2Generator(ellipse2);

    int ellipse1PointsLearn = configuration.getVariableInt("ELLIPSE_1_LEARN_POINTS");
    int ellipse2PointsLearn = configuration.getVariableInt("ELLIPSE_2_LEARN_POINTS");
    int ellipse1PointsTest = configuration.getVariableInt("ELLIPSE_1_TEST_POINTS");
    int ellipse2PointsTest = configuration.getVariableInt("ELLIPSE_2_TEST_POINTS");

    std::vector<std::vector<double>> XLearn1 = ellipse1Generator.generateLearnData(ellipse1PointsLearn, "first_ellipse");
    std::vector<std::vector<double>> XLearn2 = ellipse2Generator.generateLearnData(ellipse2PointsLearn, "second_ellipse");

//    XLearn1.push_back(std::vector<double>{0, 8});

    std::vector<double> DLearn1(XLearn1.size(), 1.0);
    std::vector<double> DLearn2(XLearn2.size(), 0.0);

    std::vector<std::vector<double>> XLearn = XLearn1 + XLearn2;
    std::vector<double> DLearn = DLearn1 + DLearn2;

    EllipseDataGenerator::writeToFile(configuration.getVariable("LEARN_DATA_FILENAME"), XLearn, DLearn);

    // Генерация данных в эллипсе для тестирования
    std::vector<std::vector<double>> XTest1 = ellipse1Generator.generateRandomData(ellipse1PointsTest);
    std::vector<std::vector<double>> XTest2 = ellipse2Generator.generateRandomData(ellipse2PointsTest);

    std::vector<double> DTest1(XTest1.size(), 1.0);
    std::vector<double> DTest2(XTest2.size(), 0.0);

    std::vector<std::vector<double>> XTest = XTest1 + XTest2;
    std::vector<double> DTest = DTest1 + DTest2;

    EllipseDataGenerator::writeToFile(configuration.getVariable("TEST_DATA_FILENAME"), XTest, DTest);

    // Создание сигмоидального нейрона, обучение, тестирование и переобучение по необходимости
    SigmoidalNeuron neuron(3, configuration);
    neuron.learn(XLearn, DLearn);

    neuron.test(XTest, DTest);
//    size_t counter_all = 0;
//    size_t counter_relearn = 0;

//    do {
//        counter_all += neuron.learn(XLearn, DLearn, 50);
//        ++counter_relearn;
//    } while (!neuron.test(XTest, DTest) && counter_relearn < C_MAX);
//
//    std::cout << "relearning cycles num = " << counter_relearn << std::endl;
//    std::cout << "all learning cycles num = " << counter_all << std::endl;

    return 0;
}
