#ifndef ELLIPSEPARAMS_H
#define ELLIPSEPARAMS_H

#include <string>

#include "config.hpp"


struct EllipseParams {
    // Полуоси эллипса
    double a;
    double b;
    // Координаты точки центра
    double x0;
    double y0;
    // Угол поворота осей эллипса
    double rotate;

    EllipseParams() = default;
    explicit EllipseParams(ConfigurationSingleton& configuration, std::string ellipseName) {
        this->a = configuration.getVariableDouble("ELLIPSE_" + ellipseName + "_A");
        this->b = configuration.getVariableDouble("ELLIPSE_" + ellipseName + "_B");
        this->x0 = configuration.getVariableDouble("ELLIPSE_" + ellipseName + "_X0");
        this->y0 = configuration.getVariableDouble("ELLIPSE_" + ellipseName + "_Y0");
        this->rotate = configuration.getVariableDouble("ELLIPSE_" + ellipseName + "_ROTATE");
    }
};

#endif // ELLIPSEPARAMS_H
