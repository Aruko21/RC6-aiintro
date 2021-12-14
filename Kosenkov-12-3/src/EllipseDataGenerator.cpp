#include <fstream>
#include <string>
#include "EllipseDataGenerator.hpp"
#include "EllipseParams.hpp"
#include <cmath>


// Сгенерерировать случайное число в интервале [l, r)
double EllipseDataGenerator::generateRandNumInRange(double l, double r) {
    double tmp = (double) std::rand() / (RAND_MAX / (r - l));
    return l + tmp;
}

std::vector<double> EllipseDataGenerator::pointToAbsoluteCoord(double xLocal, double yLocal) const {
    EllipseParams ellipse = this->ellipse;

    double x = ellipse.x0 + (xLocal * std::cos(ellipse.rotate) - yLocal * std::sin(ellipse.rotate));
    double y = ellipse.y0 + (xLocal * std::sin(ellipse.rotate) + yLocal * std::cos(ellipse.rotate));

    return std::vector<double>{x, y};
}

// Сгенерировать точки внутри эллипса
std::vector<std::vector<double>> EllipseDataGenerator::generatePointsInsideEllipse(size_t points_num) const {
    EllipseParams ellipse = this->ellipse;

    std::vector<std::vector<double>> res;

    for (size_t i = 0; i < points_num; ++i) {
        double rad = std::sqrt(generateRandNumInRange(0, 1));
        double phi = 2 * M_PI * generateRandNumInRange(0, 1);

        double x_e = rad * std::cos(phi) * ellipse.a;
        double y_e = rad * std::sin(phi) * ellipse.b;

        res.push_back(this->pointToAbsoluteCoord(x_e, y_e));
    }

    return res;
}

// Сгенерировать точки на эллипсе
std::vector<std::vector<double>> EllipseDataGenerator::generatePointsOnEllipse(size_t points_num) const {
    EllipseParams ellipse = this->ellipse;

    std::vector<std::vector<double>> result;
    double stepAngle = 360.0 / points_num;

    for (size_t i = 0; i < points_num; ++i) {
        double phi = i * stepAngle * M_PI / 180.0;

        // Переход от угла эллипса к углу образующих окружностей
        double phiAbs = std::atan2((ellipse.a * std::sin(phi)), (ellipse.b * std::cos(phi)));

        // Получение точки эллипса в локальных координатах
        double xLocal = ellipse.a * std::cos(phiAbs);
        double yLocal = ellipse.b * std::sin(phiAbs);

        result.push_back(this->pointToAbsoluteCoord(xLocal, yLocal));
    }

    return result;
}

// Генерация точек в области эллипса
std::vector<std::vector<double>> EllipseDataGenerator::generateLearnData(size_t points_num, std::string filename) const {
    // Генерация точек по периметру эллипса
    std::vector<std::vector<double>> points = generatePointsOnEllipse(points_num);

    // Генерация точек на главной оси эллипса (полуэллипс содержит главную ось)
    std::vector<double> aBegin = points.at(0); // Начало полуэллипса (начало оси)
    std::vector<double> aEnd = points.at(points_num / 2); // Конец полуэллипса (конец оси)

    // Центральная точка
    points.push_back({ (aBegin.at(0) + aEnd.at(0)) / 2, (aBegin.at(1) + aEnd.at(1)) / 2});

    auto p3 = points.at(12);

    // Точки по бокам от центра
    double x_1 = aBegin.at(0) + p3.at(0);
    double y_1 = aBegin.at(1) + p3.at(1);
    double x_2 = p3.at(0) + aEnd.at(0);
    double y_2 = p3.at(1) + aEnd.at(1);

    points.push_back({x_1 / 2, y_1 / 2});
    points.push_back({x_2 / 2, y_2 / 2});

    // Добавить несколько точек внутри эллипса
    std::vector<std::vector<double>> pointsInside = this->generatePointsInsideEllipse(1);
    points.insert(points.end(), pointsInside.begin(), pointsInside.end());

    std::string resultFilename = std::string(filename + "_" + std::to_string(this->ellipse.a) + "_" + std::to_string(this->ellipse.b) + ".dat");
    EllipseDataGenerator::writeToFile(resultFilename, points);

    return points;
}

std::vector<std::vector<double>> EllipseDataGenerator::generateRandomData(size_t points_num) const {
    return this->generatePointsInsideEllipse(points_num);
}

void EllipseDataGenerator::writeToFile(const std::string& filename, const std::vector<std::vector<double>>& res) {
    std::ofstream file(filename);
    for (auto & row : res) {
        for (auto & el : row) {
            file << el << " ";
        }
        file << std::endl;
    }
}

void EllipseDataGenerator::writeToFile(const std::string& filename,
                                  const std::vector<std::vector<double>>& res,
                                  const std::vector<double>& D_res) {
    std::ofstream file(filename);

    for (size_t i = 0; i < res.size(); i++) {
        for (auto & el : res.at(i)) {
            file << el << " ";
        }
        file << D_res.at(i) << std::endl;
    }
}

//void EllipseDataGenerator::plotDots(
//    std::vector<std::vector<double>>& XLearn,
//    std::vector<double>& DLearn,
//    std::string dataFilename,
//    std::string gnuplotFilename
//) {
//    std::ofstream gnuscript(gnuplotFilename);
//
//    gnuscript << "set term wxt title 'V(t) plot'" << std::endl;
//    gnuscript << "set xlabel 'Time t'" << std::endl;
//    gnuscript << "set ylabel 'Velocity V(t)'" << std::endl;
//    gnuscript << "set key right bottom" << std::endl;
//    gnuscript << "set grid" << std::endl;
//
//    gnuscript << "plot '" << dataFilename << "' with l lc rgb 'red' lt 1 lw 1.5 title 'Analytical Solution', \\" << std::endl;
//    gnuscript << "'" << eulerFile << "' using 1:2 with l lc rgb 'green' lt 2 lw 1.5 title 'Euler Solution', \\" << std::endl;
//    gnuscript << "'" << rungeKuttaFile << "' using 1:2 with l lc rgb 'blue' lt 4 lw 1.5 title 'Runge Kutta Solution'" << std::endl;
//    gnuscript << "pause -1" << std::endl;
//
//    system((std::string("gnuplot ") + gnuscriptFile).c_str());
//
//    for (size_t i = 0; i < XLearn.size(); i++) {
//        int color = 0;
//        int pt = 1;
//        if (DLearn.at(i) == 1) {
//            color = 2;
//            pt = 3;
//        } else {
//            color = 1;
//            pt = 2;
//        }
//        fprintf(gnuplotPipe, "plot '-' with points lc %d pt %d lw 1\n", color, pt);
//        fprintf(gnuplotPipe, "%lf %lf\n", XLearn.at(i).at(0), XLearn.at(i).at(1));
//        fprintf(gnuplotPipe, "e\n");
//    }
//    fflush(gnuplotPipe);
//}
