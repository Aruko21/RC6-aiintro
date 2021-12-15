#ifndef LAB1_INCLUDE_CONFIG_HPP
#define LAB1_INCLUDE_CONFIG_HPP

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <algorithm>

static constexpr char CONFIG_PATH[] = "config.cfg";


class ConfigurationSingleton {
private:
    std::map<std::string, std::string> config;

    ConfigurationSingleton();
    ConfigurationSingleton( const ConfigurationSingleton&);
public:
    ConfigurationSingleton& operator=(const ConfigurationSingleton& extConfig) = default;

    static ConfigurationSingleton& getInstance() {
        static ConfigurationSingleton instance;
        return instance;
    }

    std::string getVariable(std::string key);
    double getVariableDouble(std::string key);
    int getVariableInt(std::string key);
};

#endif // LAB1_INCLUDE_CONFIG_HPP