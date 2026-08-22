#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>

int main(int argc, char** argv) {
    int exit_code = 0;
    int sleep_ms = 0;
    bool split_line = false;
    std::cout << "arguments:";
    for (int i = 1; i < argc; ++i) {
        const std::string argument(argv[i]);
        std::cout << " [" << argument << "]";
        if (argument == "--exit" && i + 1 < argc) exit_code = std::stoi(argv[++i]);
        else if (argument == "--sleep" && i + 1 < argc) sleep_ms = std::stoi(argv[++i]);
        else if (argument == "--split-line") split_line = true;
    }
    std::cout << '\n';
    if (const char* value = std::getenv("SIPROSWF_TEST_VALUE")) {
        std::cout << "environment:" << value << '\n';
    }
    std::cerr << "fake-tool-stderr\n";
    std::cout.flush(); std::cerr.flush();
    if (split_line) {
        std::cout << std::string(66, '-');
        std::cout.flush();
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        std::cout << std::string(22, '-') << '\n';
        std::cout.flush();
    }
    if (sleep_ms > 0) std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
    return exit_code;
}
