#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#include "EllipseParams.hpp"


class EllipseDataGenerator {
private:
    EllipseParams ellipse;
    // Генерация точек внутри эллипса
    std::vector<std::vector<double>> generatePointsInsideEllipse(size_t points_num) const;
    // Генерация точек по периметру эллипса
    std::vector<std::vector<double>> generatePointsOnEllipse(size_t points_num) const;
    // Получение точки эллипса в абсолютных координатах (сдвиг + поворот осей)
    std::vector<double> pointToAbsoluteCoord(double xLocal, double yLocal) const;
public:
    explicit EllipseDataGenerator(EllipseParams& _ellipse) : ellipse(_ellipse) {
        srand(time(nullptr));
        this->ellipse.rotate = this->ellipse.rotate * M_PI / 180.0;
    };

    // Сгенерерировать случайное число в интервале [l, r)
    static double generateRandNumInRange(double l, double r);

    static void writeToFile(const std::string& filename, const std::vector<std::vector<double>>& res);
    static void writeToFile(const std::string& filename, const std::vector<std::vector<double>>& res, const std::vector<double>& D_res);
    static void plotDots(
        std::vector<std::vector<double>>& data,
        std::vector<double>& DLearn,
        std::string dataFilename,
        std::string gnuplotFilename
    );

    // Сгенерировать точки в эллипсе для обучения нейрона
    std::vector<std::vector<double>> generateLearnData(size_t points_num, std::string filename) const;
    std::vector<std::vector<double>> generateRandomData(size_t points_num) const;

    ~EllipseDataGenerator() = default;
};

#endif // GENERATOR_H
