#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include "Ticker.h"
#include "OptionData.h"

class FunctionsUtil
{
public:
    std::unordered_map<std::string, std::unique_ptr<Ticker>> read_csv_into_ticker_object(const std::string &fileName);
    void perform_operations(const std::unordered_map<std::string, std::unique_ptr<Ticker>> &tickers_data);
};
