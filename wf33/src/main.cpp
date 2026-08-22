#include "workflow.hpp"

#include <atomic>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#ifdef _WIN32
#  define NOMINMAX
#  include <cstdio>
#  include <windows.h>
#  include <shellapi.h>
#else
#  include <unistd.h>
#endif

namespace siproswf {

std::filesystem::path executable_path(const std::filesystem::path& fallback) {
#ifdef _WIN32
    std::wstring buffer(32768, L'\0');
    const DWORD size = GetModuleFileNameW(nullptr, buffer.data(),
                                          static_cast<DWORD>(buffer.size()));
    if (size > 0 && size < buffer.size()) {
        buffer.resize(size);
        return std::filesystem::path(buffer);
    }
#else
    std::error_code error;
    const auto path = std::filesystem::read_symlink("/proc/self/exe", error);
    if (!error) return path;
#endif
    return std::filesystem::absolute(fallback);
}

#ifdef _WIN32
static std::string utf8(const wchar_t* value) {
    if (value == nullptr || *value == L'\0') return {};
    const int size = WideCharToMultiByte(CP_UTF8, 0, value, -1, nullptr, 0,
                                         nullptr, nullptr);
    if (size <= 1) return {};
    std::string result(static_cast<std::size_t>(size), '\0');
    WideCharToMultiByte(CP_UTF8, 0, value, -1, result.data(), size,
                        nullptr, nullptr);
    result.pop_back();
    return result;
}

static std::wstring utf16(const std::string& value) {
    if (value.empty()) return {};
    const int size = MultiByteToWideChar(CP_UTF8, 0, value.c_str(),
                                         static_cast<int>(value.size()),
                                         nullptr, 0);
    if (size <= 0) return {};
    std::wstring result(static_cast<std::size_t>(size), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, value.c_str(),
                        static_cast<int>(value.size()), result.data(), size);
    return result;
}

void attach_parent_console() {
    if (!AttachConsole(ATTACH_PARENT_PROCESS)) return;
    FILE* stream = nullptr;
    freopen_s(&stream, "CONOUT$", "w", stdout);
    freopen_s(&stream, "CONOUT$", "w", stderr);
    freopen_s(&stream, "CONIN$", "r", stdin);
    SetConsoleOutputCP(CP_UTF8);
}

void configure_packaged_library_path(const std::filesystem::path& self) {
    std::error_code error;
    const std::filesystem::path library_directory = self.parent_path() / L"lib";
    if (std::filesystem::is_directory(library_directory, error)) {
        SetDllDirectoryW(library_directory.c_str());
    }
}
#endif

int run_application(const std::vector<std::string>& arguments,
                    const std::filesystem::path& self,
                    bool graphical_invocation) {
    try {
        const siproswf::ParseResult parsed = siproswf::parse_arguments(arguments);
        if (parsed.show_help) {
            siproswf::print_help(std::cout);
            return 0;
        }
        if (parsed.launch_gui) return siproswf::run_gui(self);

        std::error_code error;
        const bool existed = std::filesystem::exists(parsed.options.output, error);
        std::filesystem::create_directories(parsed.options.output);
        const std::filesystem::path log_path =
            std::filesystem::path(parsed.options.output) /
            siproswf::kWorkflowLogFilename;
        siproswf::Logger logger(log_path);
        logger.info("Workflow log: " + log_path.string());
        if (existed) logger.warning(parsed.options.output +
                                    " exists and will be overwritten");
        for (const auto& warning : parsed.warnings) logger.warning(warning);
        std::atomic_bool cancelled{false};
        try {
            siproswf::Workflow workflow(parsed.options, siproswf::locate_tools(self),
                                        logger, cancelled);
            workflow.run();
        } catch (const std::exception& error) {
            logger.error(std::string("Workflow failed: ") + error.what());
            throw;
        }
        return 0;
    } catch (const std::exception& error) {
#ifdef _WIN32
        if (graphical_invocation) {
            const std::wstring message = utf16(std::string("siproswf: ") + error.what());
            MessageBoxW(nullptr, message.c_str(), L"Sipros Workflow",
                        MB_OK | MB_ICONERROR);
        } else {
            std::cerr << "siproswf: " << error.what() << '\n';
        }
#else
        (void)graphical_invocation;
        std::cerr << "siproswf: " << error.what() << '\n';
#endif
        return 2;
    }
}

} // namespace siproswf

#ifdef _WIN32
int WINAPI wWinMain(HINSTANCE, HINSTANCE, PWSTR, int) {
    int argument_count = 0;
    LPWSTR* wide_arguments = CommandLineToArgvW(GetCommandLineW(), &argument_count);
    if (wide_arguments == nullptr || argument_count < 1) return 2;

    std::vector<std::string> arguments;
    arguments.reserve(static_cast<std::size_t>(argument_count - 1));
    for (int index = 1; index < argument_count; ++index) {
        arguments.push_back(siproswf::utf8(wide_arguments[index]));
    }
    const std::filesystem::path fallback(wide_arguments[0]);
    LocalFree(wide_arguments);

    const bool graphical_invocation = arguments.empty() ||
        (arguments.size() == 1 && arguments.front() == "--gui");
    if (!graphical_invocation) siproswf::attach_parent_console();
    const std::filesystem::path self = siproswf::executable_path(fallback);
    siproswf::configure_packaged_library_path(self);
    return siproswf::run_application(arguments, self,
                                     graphical_invocation);
}
#else
int main(int argc, char** argv) {
    std::vector<std::string> arguments;
    for (int index = 1; index < argc; ++index) arguments.emplace_back(argv[index]);
    return siproswf::run_application(
        arguments, siproswf::executable_path(argv[0]), arguments.empty());
}
#endif
