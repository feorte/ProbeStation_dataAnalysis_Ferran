#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TTree.h"

namespace {

constexpr int kHeaderLines = 30;
constexpr int kChannelHeaderLine = 15;
constexpr int kDefaultChannelCount = 264;
constexpr int kIvDummyColumns = 8;

std::string trim(const std::string &text) {
    const std::size_t first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const std::size_t last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

std::string value_after_colon(const std::string &line) {
    const std::size_t position = line.find(':');
    return position == std::string::npos ? "" : trim(line.substr(position + 1));
}

bool parse_header(std::ifstream &input, std::string &measurement_type, std::vector<int> &channels) {
    std::string line;
    for (int i = 0; i < kHeaderLines; ++i) {
        if (!std::getline(input, line)) {
            return false;
        }

        if (i == 0) {
            measurement_type = value_after_colon(line);
        } else if (i == kChannelHeaderLine) {
            const std::string channels_text = value_after_colon(line);
            if (channels_text == "None") {
                for (int channel = 1; channel <= kDefaultChannelCount; ++channel) {
                    channels.push_back(channel);
                }
            } else {
                std::stringstream stream(channels_text);
                std::string token;
                while (std::getline(stream, token, ',')) {
                    channels.push_back(std::stoi(trim(token)));
                }
            }
        }
    }
    return true;
}

struct MeasurementField {
    std::string name;
    float value = 0.0f;
    std::vector<float> values;

    MeasurementField(std::string field_name, std::size_t channel_count)
        : name(std::move(field_name)), values(channel_count, 0.0f) {}

    void add_raw_branch(TTree &tree) {
        tree.Branch(name.c_str(), &value);
    }
    void add_analysis_branch(TTree &tree) {
        tree.Branch(name.c_str(), &values);
    }
};

MeasurementField &find_field(std::vector<MeasurementField> &fields, const std::string &name) {
    for (auto &field : fields) {
        if (field.name == name) {
            return field;
        }
    }
    throw std::runtime_error("Measurement field not found: " + name);
}

bool read_measurement(std::istream &input, float &voltage, int &channel,
                      std::vector<MeasurementField> &fields, int dummy_columns) {
    if (!(input >> voltage >> channel)) {
        return false;
    }

    for (auto &field : fields) {
        if (!(input >> field.value)) {
            return false;
        }
    }

    float dummy = 0.0f;
    for (int i = 0; i < dummy_columns; ++i) {
        if (!(input >> dummy)) {
            return false;
        }
    }
    return true;
}

void compute_statistics(const std::vector<float> &values, float &mean, float &standard_deviation) {
    const float n = static_cast<float>(values.size());
    mean = 0.0f;
    for (const float value : values) {
        mean += value;
    }
    mean /= n;

    standard_deviation = 0.0f;
    for (const float value : values) {
        const float difference = value - mean;
        standard_deviation += difference * difference;
    }
    standard_deviation = std::sqrt(standard_deviation / n);
}

void process_blocks(std::istream &input, TTree &raw_tree, TTree &analysis_tree, float &voltage,
                    int &channel, std::size_t channel_count, std::vector<MeasurementField> &fields,
                    MeasurementField &temperature, MeasurementField &humidity,
                    float &mean_temperature, float &mean_humidity, float &std_temperature,
                    float &std_humidity, int dummy_columns) {
    while (true) {
        bool complete_block = true;

        for (std::size_t i = 0; i < channel_count; ++i) {
            if (!read_measurement(input, voltage, channel, fields, dummy_columns)) {
                complete_block = false;
                break;
            }

            raw_tree.Fill();
            for (auto &field : fields) {
                field.values[i] = field.value;
            }
        }

        if (!complete_block) {
            break;
        }

        compute_statistics(temperature.values, mean_temperature, std_temperature);
        compute_statistics(humidity.values, mean_humidity, std_humidity);
        analysis_tree.Fill();

        if (input.eof() || !input) {
            break;
        }
    }
}

} // namespace

int read_store(
    std::string datafile = "raw_data/CALICE_6in_256ch_CSIS2025_008_IV_CV_20250523_0_IV.txt",
    std::string sensor_number = "CSIS_008") {
    std::ifstream input(datafile);
    std::cout << "Reading file: " << datafile << std::endl;

    if (!input.is_open()) {
        std::cerr << "Could not open file: " << datafile << std::endl;
        return 1;
    }

    std::string measurement_type;
    std::vector<int> channels;
    if (!parse_header(input, measurement_type, channels)) {
        std::cerr << "Could not read the input header." << std::endl;
        return 1;
    }

    std::cout << "Measurement type: " << measurement_type << std::endl;
    std::cout << "Custom channels (" << channels.size() << "): ";
    for (const int channel : channels) {
        std::cout << channel << " ";
    }
    std::cout << std::endl;

    if (channels.empty()) {
        std::cerr << "No channels were found in the input header." << std::endl;
        return 1;
    }

    const std::size_t channel_count = channels.size();

    std::filesystem::create_directories("stored_data");
    const std::string storing_file =
        "stored_data/stored_data_" + sensor_number + "_" + measurement_type + ".root";

    auto output_file = std::make_unique<TFile>(storing_file.c_str(), "RECREATE");
    if (!output_file || output_file->IsZombie()) {
        std::cerr << "Could not create output file: " << storing_file << std::endl;
        return 1;
    }

    auto raw_tree = std::make_unique<TTree>("raw_measurements", "Raw measurements");
    auto analysis_tree = std::make_unique<TTree>("analysis", "Analysis");

    float voltage = 0.0f;
    int channel = 0;
    raw_tree->Branch("voltage", &voltage);
    raw_tree->Branch("channel", &channel);
    analysis_tree->Branch("voltage", &voltage);
    analysis_tree->Branch("channel", &channels);

    std::vector<MeasurementField> fields;
    int dummy_columns = 0;

    if (measurement_type == "IV") {
        std::cout << "Processing IV data...\n";
        fields.emplace_back("current", channel_count);
        fields.emplace_back("current_err", channel_count);
        fields.emplace_back("tot_curr", channel_count);
        fields.emplace_back("act_vlt", channel_count);
        fields.emplace_back("time", channel_count);
        fields.emplace_back("temp", channel_count);
        fields.emplace_back("hum", channel_count);
        dummy_columns = kIvDummyColumns;
    } else if (measurement_type == "CV") {
        std::cout << "Processing CV data...\n";
        fields.emplace_back("cs", channel_count);
        fields.emplace_back("cs_err", channel_count);
        fields.emplace_back("tot_curr", channel_count);
        fields.emplace_back("act_vlt", channel_count);
        fields.emplace_back("time", channel_count);
        fields.emplace_back("temp", channel_count);
        fields.emplace_back("hum", channel_count);
        fields.emplace_back("cp", channel_count);
        fields.emplace_back("cp_err", channel_count);
        fields.emplace_back("impedance", channel_count);
        fields.emplace_back("impedance_err", channel_count);
        fields.emplace_back("phase", channel_count);
        fields.emplace_back("phase_err", channel_count);
        fields.emplace_back("cs_uncorr", channel_count);
        fields.emplace_back("cp_uncorr", channel_count);
    } else {
        std::cerr << "Unknown measurement type: " << measurement_type << std::endl;
        return 1;
    }

    for (auto &field : fields) {
        field.add_raw_branch(*raw_tree);
        field.add_analysis_branch(*analysis_tree);
    }

    MeasurementField &temperature = find_field(fields, "temp");
    MeasurementField &humidity = find_field(fields, "hum");

    float mean_temperature = 0.0f;
    float mean_humidity = 0.0f;
    float std_temperature = 0.0f;
    float std_humidity = 0.0f;
    analysis_tree->Branch("mean_temp", &mean_temperature);
    analysis_tree->Branch("mean_hum", &mean_humidity);
    analysis_tree->Branch("std_temp", &std_temperature);
    analysis_tree->Branch("std_hum", &std_humidity);

    process_blocks(input, *raw_tree, *analysis_tree, voltage, channel, channel_count, fields,
                   temperature, humidity, mean_temperature, mean_humidity, std_temperature,
                   std_humidity, dummy_columns);

    raw_tree->Write();
    analysis_tree->Write();
    raw_tree.reset();
    analysis_tree.reset();
    output_file->Close();

    std::cout << "Storing file created successfully: " << storing_file << std::endl;
    return 0;
}
