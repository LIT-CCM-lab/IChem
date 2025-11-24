#ifndef IFP_MODULE_H
#define IFP_MODULE_H

#include <map>
#include <string>
#include <vector>

namespace IFPModule {

using OptionMap = std::map<std::string, std::vector<std::string>>;

void runIFP(const std::vector<std::string>& inputValues, const OptionMap& optValues);

}

#endif // IFP_MODULE_H
