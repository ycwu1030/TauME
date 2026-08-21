#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>

#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>

namespace {

std::string shell_quote(const std::string& value) {
    std::string quoted{"'"};
    for (const char character : value) {
        if (character == '\'') {
            quoted += "'\\''";
        } else {
            quoted += character;
        }
    }
    return quoted + "'";
}

std::string command(const std::string& caller_path, const std::string& input_path, const std::string& output_path,
                    const std::string& range) {
    return shell_quote(caller_path) + " --input " + shell_quote(input_path) + " --output " + shell_quote(output_path) +
           " " + range +
           " --pion-source generated --electron-polarization 0 --positron-polarization 0 --tau-mass 1.77686"
           " --electric-charge 0.313 --weak-mixing-sine 0.48 --z-mass 91.1876 --z-width 2.4952"
           " --input-sha256 a34d645cb9df4dc9a9527263277ebf3b14eb86aa9e4d74feab60744d28e84daa"
           " --implementation-revision pion-pair-delphes-export-test";
}

std::string execute(const std::string& caller_path, const std::string& input_path, const std::string& output_path,
                    const std::string& range, const std::string& standard_output_path,
                    const std::string& standard_error_path) {
    return command(caller_path, input_path, output_path, range) + " > " + shell_quote(standard_output_path) +
           " 2> " + shell_quote(standard_error_path);
}

std::string read_file(const std::string& path) {
    std::ifstream input(path);
    std::ostringstream content;
    content << input.rdbuf();
    return content.str();
}

}  // namespace

int main(int argc, char* argv[]) {
    assert(argc == 3);
    const std::string caller_path = argv[1];
    const std::string input_path = argv[2];
    const std::string output_path = "/private/tmp/tauamp_pion_pair_delphes_export_test.root";
    const std::string invalid_output_path = "/private/tmp/tauamp_pion_pair_delphes_export_invalid_range.root";
    const std::string missing_option_output_path = "/private/tmp/tauamp_pion_pair_delphes_export_missing_option.root";
    const std::string standard_output_path = "/private/tmp/tauamp_pion_pair_delphes_export.stdout";
    const std::string standard_error_path = "/private/tmp/tauamp_pion_pair_delphes_export.stderr";
    std::remove(output_path.c_str());
    std::remove(invalid_output_path.c_str());
    std::remove(missing_option_output_path.c_str());
    std::remove(standard_output_path.c_str());
    std::remove(standard_error_path.c_str());

    const int successful_run = std::system(
        execute(caller_path, input_path, output_path, "--begin 0 --end 1", standard_output_path, standard_error_path).c_str());
    assert(successful_run == 0);
    assert(read_file(standard_output_path) ==
           "entries_visited=1\n"
           "supported_input=1\n"
           "unsupported_input=0\n"
           "kinematic_no_solution=0\n"
           "kinematic_non_unique=0\n"
           "evaluable_one_solution=0\n"
           "evaluable_two_solutions=1\n");

    std::unique_ptr<TFile> output(TFile::Open(output_path.c_str(), "READ"));
    assert(output && !output->IsZombie());
    auto* events = dynamic_cast<TTree*>(output->Get("Events"));
    auto* hypotheses = dynamic_cast<TTree*>(output->Get("Hypotheses"));
    auto* manifest = dynamic_cast<TObjString*>(output->Get("TauAmpPionPairMetadata"));
    assert(events != nullptr);
    assert(hypotheses != nullptr);
    assert(manifest != nullptr);
    assert(events->GetEntries() == 1);
    assert(hypotheses->GetEntries() >= 0);
    assert(std::string(manifest->GetString().Data()).find("implementation_revision=pion-pair-delphes-export-test") !=
           std::string::npos);
    output.reset();

    const int invalid_range_run = std::system(
        execute(caller_path, input_path, invalid_output_path, "--begin 1 --end 0", standard_output_path, standard_error_path)
            .c_str());
    assert(invalid_range_run != 0);
    assert(!std::filesystem::exists(invalid_output_path));

    const std::string missing_option_command =
        shell_quote(caller_path) + " --input " + shell_quote(input_path) + " --output " +
        shell_quote(missing_option_output_path) +
        " --begin 0 --end 1 --pion-source generated --electron-polarization 0 --positron-polarization 0"
        " --tau-mass 1.77686 --electric-charge 0.313 --weak-mixing-sine 0.48 --z-mass 91.1876"
        " --input-sha256 a34d645cb9df4dc9a9527263277ebf3b14eb86aa9e4d74feab60744d28e84daa"
        " --implementation-revision pion-pair-delphes-export-test > " + shell_quote(standard_output_path) +
        " 2> " + shell_quote(standard_error_path);
    const int missing_option_run = std::system(missing_option_command.c_str());
    assert(missing_option_run != 0);
    assert(!std::filesystem::exists(missing_option_output_path));

    const int overwrite_run = std::system(
        execute(caller_path, input_path, output_path, "--begin 0 --end 1", standard_output_path, standard_error_path).c_str());
    assert(overwrite_run != 0);
    std::unique_ptr<TFile> preserved_output(TFile::Open(output_path.c_str(), "READ"));
    assert(preserved_output && !preserved_output->IsZombie());
    auto* preserved_events = dynamic_cast<TTree*>(preserved_output->Get("Events"));
    assert(preserved_events != nullptr);
    assert(preserved_events->GetEntries() == 1);

    std::remove(output_path.c_str());
    std::remove(standard_output_path.c_str());
    std::remove(standard_error_path.c_str());
}
