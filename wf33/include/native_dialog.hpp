#pragma once

#include <optional>
#include <string>

struct GLFWwindow;

namespace siproswf {

enum class NativeDialogKind {
    InputFiles,
    FastaFile,
    Directory,
};

bool native_dialog_available();
std::optional<std::string> show_native_dialog(
    GLFWwindow* owner, NativeDialogKind kind, const std::string& title);

} // namespace siproswf
