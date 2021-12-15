#include "config.hpp"

ConfigurationSingleton::ConfigurationSingleton() {
    std::ifstream file(CONFIG_PATH);

    if (!file.is_open()) {
        throw std::runtime_error("couldn't open file");
    }

    std::string line;
    while(getline(file, line)){
        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());

        if (line[0] == '#' || line.empty()) {
            continue;
        }

        size_t delimiterPos = line.find('=');
        std::string key = line.substr(0, delimiterPos);
        std::string value = line.substr(delimiterPos + 1);

        config.insert(std::pair<std::string, std::string>(key, value));
    }
}

std::string ConfigurationSingleton::getVariable(std::string key) {
    if (this->config.find(key) == this->config.end()) {
        return "";
    }

    return config[key];
}

double ConfigurationSingleton::getVariableDouble(std::string key) {
    std::string stringVariable = this->getVariable(key);
    return atof(stringVariable.c_str());
}

int ConfigurationSingleton::getVariableInt(std::string key) {
    std::string stringVariable = this->getVariable(key);
    return atoi(stringVariable.c_str());
}
