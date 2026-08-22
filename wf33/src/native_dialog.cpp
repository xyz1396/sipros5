#include "native_dialog.hpp"

#include <iterator>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef _WIN32
#  define NOMINMAX
#  define GLFW_EXPOSE_NATIVE_WIN32
#  include <windows.h>
#  include <shobjidl.h>
#  include <wrl/client.h>
#  include <GLFW/glfw3.h>
#  include <GLFW/glfw3native.h>
#endif

namespace siproswf {

#ifdef _WIN32
using Microsoft::WRL::ComPtr;

class ComApartment {
public:
    ComApartment() {
        result_ = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED | COINIT_DISABLE_OLE1DDE);
        if (FAILED(result_) && result_ != RPC_E_CHANGED_MODE) {
            throw std::runtime_error("Unable to initialize the Windows file dialog");
        }
    }

    ~ComApartment() {
        if (result_ == S_OK || result_ == S_FALSE) CoUninitialize();
    }

private:
    HRESULT result_ = E_FAIL;
};

static std::string dialog_utf8(const wchar_t* value) {
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

static std::wstring dialog_utf16(const std::string& value) {
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

std::string hresult_message(const char* operation, HRESULT result) {
    std::ostringstream message;
    message << operation << " failed (HRESULT 0x" << std::hex
            << static_cast<unsigned long>(result) << ')';
    return message.str();
}

std::optional<std::string> item_path(IShellItem* item) {
    PWSTR raw_path = nullptr;
    const HRESULT result = item->GetDisplayName(SIGDN_FILESYSPATH, &raw_path);
    if (FAILED(result)) return std::nullopt;
    const std::string path = dialog_utf8(raw_path);
    CoTaskMemFree(raw_path);
    return path;
}

#endif

bool native_dialog_available() {
#ifdef _WIN32
    return true;
#else
    return false;
#endif
}

std::optional<std::string> show_native_dialog(
    GLFWwindow* owner, NativeDialogKind kind, const std::string& title) {
#ifdef _WIN32
    ComApartment apartment;
    ComPtr<IFileOpenDialog> dialog;
    HRESULT result = CoCreateInstance(CLSID_FileOpenDialog, nullptr,
                                      CLSCTX_INPROC_SERVER,
                                      IID_PPV_ARGS(dialog.GetAddressOf()));
    if (FAILED(result)) throw std::runtime_error(hresult_message("Creating file dialog", result));

    DWORD options = 0;
    dialog->GetOptions(&options);
    options |= FOS_FORCEFILESYSTEM | FOS_PATHMUSTEXIST | FOS_NOCHANGEDIR |
               FOS_DONTADDTORECENT;
    if (kind == NativeDialogKind::Directory) options |= FOS_PICKFOLDERS;
    else options |= FOS_FILEMUSTEXIST;
    if (kind == NativeDialogKind::InputFiles) options |= FOS_ALLOWMULTISELECT;
    dialog->SetOptions(options);

    if (kind == NativeDialogKind::InputFiles) {
        static const COMDLG_FILTERSPEC filters[] = {
            {L"Mass spectrometry inputs", L"*.raw;*.h5;*.hdf5;*.zip"},
            {L"All files", L"*.*"},
        };
        dialog->SetFileTypes(static_cast<UINT>(std::size(filters)), filters);
        dialog->SetFileTypeIndex(1);
    } else if (kind == NativeDialogKind::FastaFile) {
        static const COMDLG_FILTERSPEC filters[] = {
            {L"FASTA databases", L"*.faa;*.fasta;*.fa;*.fas"},
            {L"All files", L"*.*"},
        };
        dialog->SetFileTypes(static_cast<UINT>(std::size(filters)), filters);
        dialog->SetFileTypeIndex(1);
    }

    const std::wstring requested_title = dialog_utf16(title);
    if (!requested_title.empty()) dialog->SetTitle(requested_title.c_str());
    else if (kind == NativeDialogKind::InputFiles) dialog->SetTitle(L"Select raw or HDF5 input files");
    else if (kind == NativeDialogKind::FastaFile) dialog->SetTitle(L"Select a FASTA database");
    else dialog->SetTitle(L"Select a directory");

    result = dialog->Show(glfwGetWin32Window(owner));
    if (result == HRESULT_FROM_WIN32(ERROR_CANCELLED)) return std::nullopt;
    if (FAILED(result)) throw std::runtime_error(hresult_message("Opening file dialog", result));

    if (kind == NativeDialogKind::InputFiles) {
        ComPtr<IShellItemArray> items;
        result = dialog->GetResults(items.GetAddressOf());
        if (FAILED(result)) throw std::runtime_error(hresult_message("Reading file selection", result));
        DWORD count = 0;
        items->GetCount(&count);
        std::ostringstream joined;
        bool first = true;
        for (DWORD index = 0; index < count; ++index) {
            ComPtr<IShellItem> item;
            if (FAILED(items->GetItemAt(index, item.GetAddressOf()))) continue;
            const auto path = item_path(item.Get());
            if (!path || path->empty()) continue;
            if (!first) joined << ',';
            joined << *path;
            first = false;
        }
        if (first) return std::nullopt;
        return joined.str();
    }

    ComPtr<IShellItem> item;
    result = dialog->GetResult(item.GetAddressOf());
    if (FAILED(result)) throw std::runtime_error(hresult_message("Reading path selection", result));
    return item_path(item.Get());
#else
    (void)owner;
    (void)kind;
    (void)title;
    return std::nullopt;
#endif
}

} // namespace siproswf
