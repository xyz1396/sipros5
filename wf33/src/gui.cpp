#include "workflow.hpp"

#include "native_dialog.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <filesystem>
#include <mutex>
#include <optional>
#include <sstream>
#include <string>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>

#ifdef _WIN32
#  define NOMINMAX
#  define GLFW_EXPOSE_NATIVE_WIN32
#  include <windows.h>
#endif
#include <GLFW/glfw3.h>
#ifdef _WIN32
#  include <GLFW/glfw3native.h>
#endif
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <imgui_stdlib.h>

namespace siproswf {

constexpr float kBaseFontSize = 20.0f;
constexpr float kBaseTitleHeight = 52.0f;
constexpr float kBaseTitleButtonWidth = 46.0f;
constexpr float kBaseResizeBorder = 7.0f;
constexpr int kBaseMinimumWidth = 760;
constexpr int kBaseMinimumHeight = 560;
constexpr int kSiprosIconResource = 101;
constexpr int kRegularDefaultMaxPtmCount = 3;
constexpr float kMaximumDpiScale = 5.0f;
constexpr float kDpiScaleEpsilon = 0.01f;
constexpr std::size_t kMaximumGuiLogLines = 10000;
constexpr double kActiveEventWaitSeconds = 0.50;

struct PtmChoice {
    const char* selector;
    const char* label;
};

constexpr std::array<PtmChoice, 13> kVariablePtmChoices{{
    {"oxidation", "Oxidation (M)"},
    {"deamidation", "Deamidation (N/Q)"},
    {"phosphorylation", "Phosphorylation (S/T/Y/H/D)"},
    {"phosphorylation-loss-hpo3", "Phosphorylation, HPO3 loss"},
    {"phosphorylation-loss-hpo3-h2o", "Phosphorylation, HPO3 + H2O loss"},
    {"acetylation", "Acetylation (K)"},
    {"mono-methylation", "Mono-methylation (K/R/E/D)"},
    {"di-methylation", "Di-methylation (K/R)"},
    {"tri-methylation", "Tri-methylation (K)"},
    {"s-nitrosylation", "S-nitrosylation (C)"},
    {"nitration", "Nitration (Y)"},
    {"iaa-blocking", "IAA blocking (C)"},
    {"beta-methylthiolation", "Beta-methylthiolation (D)"},
}};

enum ResizeEdge {
    ResizeNone = 0,
    ResizeLeft = 1 << 0,
    ResizeRight = 1 << 1,
    ResizeTop = 1 << 2,
    ResizeBottom = 1 << 3,
};

enum class PathTarget {
    InputFiles,
    InputDirectory,
    FastaFile,
    OutputDirectory,
    NegativeControlFiles,
};

struct GuiState {
    std::mutex mutex;
    std::deque<std::string> logs;
    std::size_t hidden_log_lines = 0;
    std::string status = "Ready";
    std::thread worker;
    std::atomic_bool running{false};
    std::atomic_bool cancelled{false};
    bool scroll_to_bottom = false;
};

struct FramelessState {
    bool previous_left_pressed = false;
    bool dragging = false;
    int resize_edges = ResizeNone;
    int start_x = 0;
    int start_y = 0;
    int start_width = 0;
    int start_height = 0;
    double start_global_x = 0.0;
    double start_global_y = 0.0;
    double drag_offset_x = 0.0;
    double drag_offset_y = 0.0;
};

struct BrowserState {
    bool request_open = false;
    bool show_all_files = false;
    PathTarget target = PathTarget::InputFiles;
    std::filesystem::path current_directory;
    std::filesystem::path selected_path;
    std::string location;
    std::string error;
};

struct BrowserEntry {
    std::filesystem::path path;
    bool directory = false;
    std::uintmax_t size = 0;
};

struct PtmSelection {
    std::array<bool, kVariablePtmChoices.size()> variable{};
    bool fixed_carbamidomethyl = true;
};

struct ExitState {
    bool confirmation_requested = false;
    bool confirmed = false;
};

struct WindowState {
    ExitState exit;
    float content_scale = 1.0f;
    bool dpi_change_pending = false;
    GLFWwindowcontentscalefun previous_content_scale_callback = nullptr;
};

struct MonitorWorkArea {
    int x = 0;
    int y = 0;
    int width = 1024;
    int height = 768;
};

enum class TitleControl {
    Minimize,
    Maximize,
    Close,
};

bool initialize_glfw(bool& use_custom_frame) {
#if defined(__linux__)
    use_custom_frame = false;
#  if defined(GLFW_PLATFORM_X11)
    // GLFW otherwise prefers Wayland when both WSLg display servers are
    // available.  Wayland requires EGL and does not allow glfwSetWindowPos,
    // while the X11/GLX path supports the custom movable frame used here.
    const char* display = std::getenv("DISPLAY");
    if (display != nullptr && *display != '\0' &&
        glfwPlatformSupported(GLFW_PLATFORM_X11) == GLFW_TRUE) {
        glfwInitHint(GLFW_PLATFORM, GLFW_PLATFORM_X11);
        if (glfwInit() == GLFW_TRUE) {
            use_custom_frame = true;
            return true;
        }
        glfwInitHint(GLFW_PLATFORM, GLFW_ANY_PLATFORM);
    }
#  endif
    return glfwInit() == GLFW_TRUE;
#else
    use_custom_frame = true;
    return glfwInit() == GLFW_TRUE;
#endif
}

static std::string lower_copy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

static std::string trim_copy(const std::string& value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string sample_name_from_path(const std::filesystem::path& path) {
    const std::string filename = path.filename().string();
    const std::string lowered = lower_copy(filename);
    if (lowered.size() >= 6 && lowered.rfind(".d.zip") == lowered.size() - 6) {
        return filename.substr(0, filename.size() - 6);
    }
    return path.stem().string();
}

std::string sample_names_from_paths(const std::string& selected_paths) {
    std::istringstream input(selected_paths);
    std::vector<std::string> names;
    std::string path;
    while (std::getline(input, path, ',')) {
        const std::string name = sample_name_from_path(std::filesystem::path(trim_copy(path)));
        if (!name.empty() && std::find(names.begin(), names.end(), name) == names.end()) {
            names.push_back(name);
        }
    }
    std::ostringstream joined;
    for (std::size_t index = 0; index < names.size(); ++index) {
        if (index != 0) joined << ',';
        joined << names[index];
    }
    return joined.str();
}

#ifdef _WIN32
std::vector<unsigned char> icon_pixels(int dimension) {
    HICON icon = static_cast<HICON>(LoadImageW(
        GetModuleHandleW(nullptr), MAKEINTRESOURCEW(kSiprosIconResource),
        IMAGE_ICON, dimension, dimension, LR_DEFAULTCOLOR));
    if (icon == nullptr) return {};

    BITMAPINFO bitmap_info{};
    bitmap_info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bitmap_info.bmiHeader.biWidth = dimension;
    bitmap_info.bmiHeader.biHeight = -dimension;
    bitmap_info.bmiHeader.biPlanes = 1;
    bitmap_info.bmiHeader.biBitCount = 32;
    bitmap_info.bmiHeader.biCompression = BI_RGB;

    HDC screen = GetDC(nullptr);
    HDC memory = CreateCompatibleDC(screen);
    void* raw_pixels = nullptr;
    HBITMAP bitmap = CreateDIBSection(screen, &bitmap_info, DIB_RGB_COLORS,
                                      &raw_pixels, nullptr, 0);
    std::vector<unsigned char> rgba;
    if (memory != nullptr && bitmap != nullptr && raw_pixels != nullptr) {
        const HGDIOBJ previous = SelectObject(memory, bitmap);
        const std::size_t byte_count = static_cast<std::size_t>(dimension) *
                                       static_cast<std::size_t>(dimension) * 4;
        std::memset(raw_pixels, 0, byte_count);
        if (DrawIconEx(memory, 0, 0, icon, dimension, dimension, 0, nullptr,
                       DI_NORMAL) != FALSE) {
            const auto* bgra = static_cast<const unsigned char*>(raw_pixels);
            rgba.resize(byte_count);
            for (std::size_t offset = 0; offset < byte_count; offset += 4) {
                const unsigned int alpha = bgra[offset + 3];
                auto unpremultiply = [alpha](unsigned int channel) {
                    return static_cast<unsigned char>(
                        alpha > 0 && alpha < 255
                            ? std::min(255U, (channel * 255U + alpha / 2U) / alpha)
                            : channel);
                };
                rgba[offset] = unpremultiply(bgra[offset + 2]);
                rgba[offset + 1] = unpremultiply(bgra[offset + 1]);
                rgba[offset + 2] = unpremultiply(bgra[offset]);
                rgba[offset + 3] = static_cast<unsigned char>(alpha);
            }
        }
        SelectObject(memory, previous);
    }
    if (bitmap != nullptr) DeleteObject(bitmap);
    if (memory != nullptr) DeleteDC(memory);
    if (screen != nullptr) ReleaseDC(nullptr, screen);
    DestroyIcon(icon);
    return rgba;
}

void install_native_window_icon(GLFWwindow* window) {
    HWND handle = glfwGetWin32Window(window);
    HINSTANCE instance = GetModuleHandleW(nullptr);
    HICON small_icon = static_cast<HICON>(LoadImageW(
        instance, MAKEINTRESOURCEW(kSiprosIconResource), IMAGE_ICON,
        GetSystemMetrics(SM_CXSMICON), GetSystemMetrics(SM_CYSMICON), LR_SHARED));
    HICON large_icon = static_cast<HICON>(LoadImageW(
        instance, MAKEINTRESOURCEW(kSiprosIconResource), IMAGE_ICON,
        GetSystemMetrics(SM_CXICON), GetSystemMetrics(SM_CYICON), LR_SHARED));
    if (small_icon != nullptr) SendMessageW(handle, WM_SETICON, ICON_SMALL,
                                            reinterpret_cast<LPARAM>(small_icon));
    if (large_icon != nullptr) SendMessageW(handle, WM_SETICON, ICON_BIG,
                                            reinterpret_cast<LPARAM>(large_icon));

    constexpr int icon_size = 64;
    std::vector<unsigned char> rgba = icon_pixels(icon_size);
    if (!rgba.empty()) {
        GLFWimage image{icon_size, icon_size, rgba.data()};
        glfwSetWindowIcon(window, 1, &image);
    }
}

unsigned int create_title_icon_texture(int dimension) {
    const std::vector<unsigned char> rgba = icon_pixels(dimension);
    if (rgba.empty()) return 0;
    GLint previous_texture = 0;
    glGetIntegerv(GL_TEXTURE_BINDING_2D, &previous_texture);
    unsigned int texture = 0;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, dimension, dimension, 0,
                 GL_RGBA, GL_UNSIGNED_BYTE, rgba.data());
    glBindTexture(GL_TEXTURE_2D, static_cast<unsigned int>(previous_texture));
    return texture;
}
#else
void install_native_window_icon(GLFWwindow*) {}
unsigned int create_title_icon_texture(int) { return 0; }
#endif

void append_log(GuiState& state, const std::string& message) {
    if (trim_copy(message).empty()) return;
    std::lock_guard<std::mutex> lock(state.mutex);
    if (state.logs.size() == kMaximumGuiLogLines) {
        state.logs.pop_front();
        ++state.hidden_log_lines;
    }
    state.logs.push_back(message);
    state.scroll_to_bottom = true;
}

void set_status(GuiState& state, const std::string& status) {
    bool changed = false;
    {
        std::lock_guard<std::mutex> lock(state.mutex);
        changed = state.status != status;
        state.status = status;
    }
    if (changed) glfwPostEmptyEvent();
}

std::string current_status(GuiState& state) {
    std::lock_guard<std::mutex> lock(state.mutex);
    return state.status;
}

Logger::Sink gui_log_sink(GuiState& state) {
    return [&state](const std::string& line) {
        append_log(state, line);
        for (const auto* marker : {
                 "Fast SIP phase", "Preparing Raxport",
                 "cache preparation", "paired target/decoy",
                 "Running Aerith"}) {
            const auto position = line.find(marker);
            if (position != std::string::npos) {
                set_status(state, line.substr(position));
                break;
            }
        }
    };
}

std::string preview(const WorkflowOptions& options) {
    std::ostringstream command;
    command << "siproswf ";
    if (options.mode == WorkflowMode::RegularFasta) command << "--regular-fasta-search ";
    else if (options.mode == WorkflowMode::SipFasta) command << "--sip-fasta-search ";
    else command << "--fast-sip-search ";
    command << "-i " << quote_display_argument(options.input)
            << " -f " << quote_display_argument(options.fasta)
            << " -o " << quote_display_argument(options.output);
    if (options.mode != WorkflowMode::RegularFasta) command << " -e " << options.element;
    if (options.sip_range) command << " -r " << quote_display_argument(*options.sip_range);
    if (options.precision) command << " -p " << quote_display_argument(*options.precision);
    command << " -t " << options.threads
            << " --toleranceMS1 " << options.tolerance_ms1
            << " --toleranceMS2 " << options.tolerance_ms2
            << " --nPrecursor " << options.n_precursor
            << " --product-top-isotopes " << options.product_top_isotopes
            << " --aerith-sample-parallelism " << options.aerith_sample_parallelism
            << " --topN " << options.top_psms_per_scan
            << " --rt-tolerance " << options.rt_tolerance
            << " --sfi-envelope-top-n " << options.sfi_envelope_top_n
            << " --mvh-cascade-top-n " << options.mvh_cascade_top_n;
    for (const auto& ptm : options.ptms) command << " --ptm " << quote_display_argument(ptm);
    for (const auto& ptm : options.fixed_ptms) command << " --fixed-ptm " << quote_display_argument(ptm);
    if (options.max_ptm_count) command << " --max-ptm-count " << *options.max_ptm_count;
    if (options.mode != WorkflowMode::RegularFasta &&
        !options.negative_control.empty()) {
        command << " --negative_control " << quote_display_argument(options.negative_control)
                << " --label_threshold " << options.label_threshold;
    }
    if (options.ignore_pct) command << " --ignorePCT";
    if (options.dry_run) command << " --dryrun";
    return command.str();
}

int resize_edges_at(double cursor_x, double cursor_y, int width, int height,
                    float border) {
    int edges = ResizeNone;
    if (cursor_x >= 0.0 && cursor_x <= border) edges |= ResizeLeft;
    else if (cursor_x >= static_cast<double>(width) - border &&
             cursor_x <= static_cast<double>(width)) edges |= ResizeRight;
    if (cursor_y >= 0.0 && cursor_y <= border) edges |= ResizeTop;
    else if (cursor_y >= static_cast<double>(height) - border &&
             cursor_y <= static_cast<double>(height)) edges |= ResizeBottom;
    return edges;
}

int update_frameless_interaction(GLFWwindow* window, FramelessState& state,
                                 float title_height, float title_controls_width,
                                 float resize_border, int minimum_width,
                                 int minimum_height) {
    int window_x = 0, window_y = 0, width = 0, height = 0;
    double cursor_x = 0.0, cursor_y = 0.0;
    glfwGetWindowPos(window, &window_x, &window_y);
    glfwGetWindowSize(window, &width, &height);
    glfwGetCursorPos(window, &cursor_x, &cursor_y);
    const bool maximized = glfwGetWindowAttrib(window, GLFW_MAXIMIZED) == GLFW_TRUE;
    const bool left_pressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    const bool just_pressed = left_pressed && !state.previous_left_pressed;
    const double global_x = static_cast<double>(window_x) + cursor_x;
    const double global_y = static_cast<double>(window_y) + cursor_y;
    const int hovered_edges = maximized
        ? ResizeNone
        : resize_edges_at(cursor_x, cursor_y, width, height, resize_border);

    if (just_pressed && !maximized) {
        if (hovered_edges != ResizeNone) {
            state.resize_edges = hovered_edges;
            state.start_x = window_x;
            state.start_y = window_y;
            state.start_width = width;
            state.start_height = height;
            state.start_global_x = global_x;
            state.start_global_y = global_y;
        } else if (cursor_y >= resize_border && cursor_y < title_height &&
                   cursor_x < static_cast<double>(width) - title_controls_width) {
            state.dragging = true;
            state.drag_offset_x = cursor_x;
            state.drag_offset_y = cursor_y;
        }
    }

    if (!left_pressed) {
        state.dragging = false;
        state.resize_edges = ResizeNone;
    } else if (state.resize_edges != ResizeNone) {
        const int delta_x = static_cast<int>(std::lround(global_x - state.start_global_x));
        const int delta_y = static_cast<int>(std::lround(global_y - state.start_global_y));
        int new_x = state.start_x;
        int new_y = state.start_y;
        int new_width = state.start_width;
        int new_height = state.start_height;
        if ((state.resize_edges & ResizeLeft) != 0) {
            new_x = state.start_x + delta_x;
            new_width = state.start_width - delta_x;
            if (new_width < minimum_width) {
                new_width = minimum_width;
                new_x = state.start_x + state.start_width - minimum_width;
            }
        } else if ((state.resize_edges & ResizeRight) != 0) {
            new_width = std::max(minimum_width, state.start_width + delta_x);
        }
        if ((state.resize_edges & ResizeTop) != 0) {
            new_y = state.start_y + delta_y;
            new_height = state.start_height - delta_y;
            if (new_height < minimum_height) {
                new_height = minimum_height;
                new_y = state.start_y + state.start_height - minimum_height;
            }
        } else if ((state.resize_edges & ResizeBottom) != 0) {
            new_height = std::max(minimum_height, state.start_height + delta_y);
        }
        if (new_x != window_x || new_y != window_y) glfwSetWindowPos(window, new_x, new_y);
        if (new_width != width || new_height != height) glfwSetWindowSize(window, new_width, new_height);
    } else if (state.dragging) {
        glfwSetWindowPos(window,
                         static_cast<int>(std::lround(global_x - state.drag_offset_x)),
                         static_cast<int>(std::lround(global_y - state.drag_offset_y)));
    }

    state.previous_left_pressed = left_pressed;
    return state.resize_edges != ResizeNone ? state.resize_edges : hovered_edges;
}

void set_resize_cursor(int edges) {
    const bool horizontal = (edges & (ResizeLeft | ResizeRight)) != 0;
    const bool vertical = (edges & (ResizeTop | ResizeBottom)) != 0;
    if (horizontal && vertical) {
        const bool northwest_southeast =
            ((edges & ResizeLeft) != 0 && (edges & ResizeTop) != 0) ||
            ((edges & ResizeRight) != 0 && (edges & ResizeBottom) != 0);
        ImGui::SetMouseCursor(northwest_southeast
                                  ? ImGuiMouseCursor_ResizeNWSE
                                  : ImGuiMouseCursor_ResizeNESW);
    } else if (horizontal) {
        ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
    } else if (vertical) {
        ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeNS);
    }
}

bool title_control(ImDrawList* draw_list, const char* id, const ImVec2& position,
                   const ImVec2& size, TitleControl control, bool maximized,
                   float scale) {
    ImGui::SetCursorScreenPos(position);
    ImGui::InvisibleButton(id, size);
    const bool hovered = ImGui::IsItemHovered();
    const bool held = ImGui::IsItemActive();
    if (hovered || held) {
        const ImU32 color = control == TitleControl::Close
            ? IM_COL32(196, 52, 58, held ? 255 : 225)
            : ImGui::GetColorU32(held ? ImGuiCol_ButtonActive : ImGuiCol_ButtonHovered);
        draw_list->AddRectFilled(position,
                                 ImVec2(position.x + size.x, position.y + size.y),
                                 color);
    }

    const ImU32 icon_color = ImGui::GetColorU32(ImGuiCol_Text);
    const float thickness = std::max(1.0f, 1.35f * scale);
    const ImVec2 center(position.x + size.x * 0.5f,
                        position.y + size.y * 0.5f);
    const float radius = 6.0f * scale;
    if (control == TitleControl::Minimize) {
        draw_list->AddLine(ImVec2(center.x - radius, center.y + 3.0f * scale),
                           ImVec2(center.x + radius, center.y + 3.0f * scale),
                           icon_color, thickness);
    } else if (control == TitleControl::Close) {
        draw_list->AddLine(ImVec2(center.x - radius, center.y - radius),
                           ImVec2(center.x + radius, center.y + radius),
                           icon_color, thickness);
        draw_list->AddLine(ImVec2(center.x + radius, center.y - radius),
                           ImVec2(center.x - radius, center.y + radius),
                           icon_color, thickness);
    } else if (maximized) {
        const float offset = 2.5f * scale;
        draw_list->AddRect(ImVec2(center.x - radius + offset, center.y - radius - offset),
                           ImVec2(center.x + radius + offset, center.y + radius - offset),
                           icon_color, 0.0f, 0, thickness);
        draw_list->AddRectFilled(ImVec2(center.x - radius - offset, center.y - radius + offset),
                                 ImVec2(center.x + radius - offset, center.y + radius + offset),
                                 ImGui::GetColorU32(ImGuiCol_TitleBg));
        draw_list->AddRect(ImVec2(center.x - radius - offset, center.y - radius + offset),
                           ImVec2(center.x + radius - offset, center.y + radius + offset),
                           icon_color, 0.0f, 0, thickness);
    } else {
        draw_list->AddRect(ImVec2(center.x - radius, center.y - radius),
                           ImVec2(center.x + radius, center.y + radius),
                           icon_color, 0.0f, 0, thickness);
    }
    return ImGui::IsItemClicked(ImGuiMouseButton_Left);
}

void render_title_bar(GLFWwindow* window, float scale, float title_height,
                      float button_width, unsigned int title_icon_texture,
                      ExitState& exit_state) {
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    const ImVec2 origin = ImGui::GetWindowPos();
    const float width = ImGui::GetWindowWidth();
    const ImVec2 bar_end(origin.x + width, origin.y + title_height);
    draw_list->AddRectFilled(origin, bar_end, ImGui::GetColorU32(ImGuiCol_TitleBg));
    draw_list->AddLine(ImVec2(origin.x, bar_end.y - 1.0f),
                       ImVec2(bar_end.x, bar_end.y - 1.0f),
                       ImGui::GetColorU32(ImGuiCol_Border));

    float text_x = origin.x + 16.0f * scale;
    if (title_icon_texture != 0) {
        const float icon_size = title_height - 16.0f * scale;
        const ImVec2 icon_min(origin.x + 10.0f * scale,
                              origin.y + (title_height - icon_size) * 0.5f);
        const ImVec2 icon_max(icon_min.x + icon_size, icon_min.y + icon_size);
        draw_list->AddImage(static_cast<ImTextureID>(title_icon_texture),
                            icon_min, icon_max);
        text_x = icon_max.x + 9.0f * scale;
    }
    ImGui::SetCursorScreenPos(ImVec2(
        text_x, origin.y + (title_height - ImGui::GetTextLineHeight()) * 0.5f));
    ImGui::PushFont(nullptr, ImGui::GetStyle().FontSizeBase * 1.08f);
    ImGui::TextUnformatted("Sipros Workflow");
    ImGui::PopFont();

    const bool maximized = glfwGetWindowAttrib(window, GLFW_MAXIMIZED) == GLFW_TRUE;
    const float button_height = title_height;
    float button_x = origin.x + width - 3.0f * button_width;
    if (title_control(draw_list, "##minimize", ImVec2(button_x, origin.y),
                      ImVec2(button_width, button_height), TitleControl::Minimize,
                      maximized, scale)) {
        glfwIconifyWindow(window);
    }
    button_x += button_width;
    if (title_control(draw_list, "##maximize", ImVec2(button_x, origin.y),
                      ImVec2(button_width, button_height), TitleControl::Maximize,
                      maximized, scale)) {
        if (maximized) glfwRestoreWindow(window);
        else glfwMaximizeWindow(window);
    }
    button_x += button_width;
    if (title_control(draw_list, "##close", ImVec2(button_x, origin.y),
                      ImVec2(button_width, button_height), TitleControl::Close,
                      maximized, scale)) {
        exit_state.confirmation_requested = true;
    }
}

float normalized_content_scale(float scale_x, float scale_y) {
    const float scale = std::max(scale_x, scale_y);
    if (!std::isfinite(scale) || scale <= 0.0f) return 1.0f;
    return std::clamp(scale, 1.0f, kMaximumDpiScale);
}

void handle_window_content_scale(GLFWwindow* window, float scale_x,
                                 float scale_y) {
    auto* state = static_cast<WindowState*>(glfwGetWindowUserPointer(window));
    if (state == nullptr) return;
    state->content_scale = normalized_content_scale(scale_x, scale_y);
    state->dpi_change_pending = true;
    if (state->previous_content_scale_callback != nullptr &&
        state->previous_content_scale_callback != handle_window_content_scale) {
        state->previous_content_scale_callback(window, scale_x, scale_y);
    }
}

void handle_window_close(GLFWwindow* window) {
    auto* state = static_cast<WindowState*>(glfwGetWindowUserPointer(window));
    if (state == nullptr || state->exit.confirmed) return;
    state->exit.confirmation_requested = true;
    glfwSetWindowShouldClose(window, GLFW_FALSE);
}

void render_exit_confirmation(GLFWwindow* window, ExitState& exit_state,
                              const GuiState& workflow_state, float scale) {
    if (exit_state.confirmation_requested) {
        ImGui::OpenPopup("Exit Sipros Workflow?##exit_confirmation");
        exit_state.confirmation_requested = false;
    }
    const ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f),
                            ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(430.0f * scale, 0.0f), ImGuiCond_Appearing);
    if (!ImGui::BeginPopupModal("Exit Sipros Workflow?##exit_confirmation", nullptr,
                                ImGuiWindowFlags_AlwaysAutoResize |
                                    ImGuiWindowFlags_NoSavedSettings)) {
        return;
    }
    ImGui::TextWrapped("Are you sure you want to exit Sipros Workflow?");
    if (workflow_state.running.load()) {
        ImGui::Spacing();
        ImGui::TextColored(ImVec4(1.0f, 0.72f, 0.30f, 1.0f),
                           "The running workflow will be stopped.");
    }
    ImGui::Spacing();
    const float button_width = 120.0f * scale;
    if (ImGui::Button("Exit", ImVec2(button_width, 0.0f))) {
        exit_state.confirmed = true;
        glfwSetWindowShouldClose(window, GLFW_TRUE);
        ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel", ImVec2(button_width, 0.0f))) {
        ImGui::CloseCurrentPopup();
    }
    ImGui::SetItemDefaultFocus();
    ImGui::EndPopup();
}

void assign_path(PathTarget target, const std::string& path, WorkflowOptions& form) {
    switch (target) {
        case PathTarget::InputFiles:
        case PathTarget::InputDirectory: form.input = path; break;
        case PathTarget::FastaFile: form.fasta = path; break;
        case PathTarget::OutputDirectory: form.output = path; break;
        case PathTarget::NegativeControlFiles:
            form.negative_control = sample_names_from_paths(path);
            break;
    }
}

std::string value_for_target(PathTarget target, const WorkflowOptions& form) {
    switch (target) {
        case PathTarget::InputFiles:
        case PathTarget::InputDirectory: return form.input;
        case PathTarget::FastaFile: return form.fasta;
        case PathTarget::OutputDirectory: return form.output;
        case PathTarget::NegativeControlFiles: return {};
    }
    return {};
}

bool browser_accepts_file(PathTarget target, const std::filesystem::path& path) {
    const std::string name = lower_copy(path.filename().string());
    if (target == PathTarget::FastaFile) {
        return name.size() >= 3 &&
               (name.rfind(".faa") == name.size() - 4 ||
                name.rfind(".fasta") == name.size() - 6 ||
                name.rfind(".fa") == name.size() - 3 ||
                name.rfind(".fas") == name.size() - 4);
    }
    if (target == PathTarget::InputFiles ||
        target == PathTarget::NegativeControlFiles) {
        return (name.size() >= 4 && name.rfind(".raw") == name.size() - 4) ||
               (name.size() >= 3 && name.rfind(".h5") == name.size() - 3) ||
               (name.size() >= 5 && name.rfind(".hdf5") == name.size() - 5) ||
               (name.size() >= 6 && name.rfind(".d.zip") == name.size() - 6);
    }
    return false;
}

void initialize_browser(BrowserState& browser, PathTarget target,
                        const WorkflowOptions& form) {
    browser.target = target;
    browser.request_open = true;
    browser.selected_path.clear();
    browser.error.clear();
    std::error_code error;
    std::filesystem::path start;
    std::string current_value = value_for_target(target, form);
    const auto comma = current_value.find(',');
    if (comma != std::string::npos) current_value.resize(comma);
    if (!current_value.empty()) {
        std::filesystem::path candidate(current_value);
        if (std::filesystem::is_directory(candidate, error)) start = candidate;
        else if (candidate.has_parent_path() && std::filesystem::is_directory(candidate.parent_path(), error)) {
            start = candidate.parent_path();
        }
    }
    if (start.empty()) start = std::filesystem::current_path(error);
    if (start.empty()) start = std::filesystem::path(".");
    browser.current_directory = std::filesystem::absolute(start, error).lexically_normal();
    browser.location = browser.current_directory.string();
}

void request_path(GLFWwindow* window, PathTarget target, WorkflowOptions& form,
                  BrowserState& browser, GuiState& state) {
    if (!native_dialog_available()) {
        initialize_browser(browser, target, form);
        return;
    }
    NativeDialogKind kind = NativeDialogKind::Directory;
    if (target == PathTarget::InputFiles ||
        target == PathTarget::NegativeControlFiles) {
        kind = NativeDialogKind::InputFiles;
    }
    else if (target == PathTarget::FastaFile) kind = NativeDialogKind::FastaFile;
    try {
        const char* title = target == PathTarget::NegativeControlFiles
            ? "Select negative-control sample files"
            : target == PathTarget::InputFiles
                ? "Select raw or HDF5 input files"
                : target == PathTarget::FastaFile
                    ? "Select a FASTA database"
                    : target == PathTarget::OutputDirectory
                        ? "Select the workflow output directory"
                        : "Select the raw-data input directory";
        const auto selected = show_native_dialog(window, kind, title);
        if (selected) assign_path(target, *selected, form);
    } catch (const std::exception& error) {
        append_log(state, std::string("File dialog error: ") + error.what());
        set_status(state, "File selection failed");
    }
}

void navigate_browser(BrowserState& browser, const std::filesystem::path& directory) {
    std::error_code error;
    const std::filesystem::path resolved = std::filesystem::absolute(directory, error).lexically_normal();
    if (error || !std::filesystem::is_directory(resolved, error)) {
        browser.error = "Directory is not accessible: " + directory.string();
        return;
    }
    browser.current_directory = resolved;
    browser.location = resolved.string();
    browser.selected_path.clear();
    browser.error.clear();
}

std::string display_size(std::uintmax_t bytes) {
    std::ostringstream output;
    if (bytes >= 1024ULL * 1024ULL * 1024ULL) {
        output.setf(std::ios::fixed); output.precision(1);
        output << static_cast<double>(bytes) / (1024.0 * 1024.0 * 1024.0) << " GB";
    } else if (bytes >= 1024ULL * 1024ULL) {
        output.setf(std::ios::fixed); output.precision(1);
        output << static_cast<double>(bytes) / (1024.0 * 1024.0) << " MB";
    } else if (bytes >= 1024ULL) {
        output << (bytes / 1024ULL) << " KB";
    } else {
        output << bytes << " B";
    }
    return output.str();
}

void render_browser(BrowserState& browser, WorkflowOptions& form, float scale) {
    if (browser.request_open) {
        ImGui::OpenPopup("Select workflow path##file_browser");
        browser.request_open = false;
    }
    const ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowSize(ImVec2(std::min(820.0f * scale, io.DisplaySize.x * 0.92f),
                                    std::min(600.0f * scale, io.DisplaySize.y * 0.90f)),
                             ImGuiCond_Appearing);
    if (!ImGui::BeginPopupModal("Select workflow path##file_browser", nullptr,
                                ImGuiWindowFlags_NoCollapse)) return;

    ImGui::SetNextItemWidth(-150.0f * scale);
    const bool enter_location = ImGui::InputText("##browser_location", &browser.location,
                                                 ImGuiInputTextFlags_EnterReturnsTrue);
    ImGui::SameLine();
    if (ImGui::Button("Go", ImVec2(54.0f * scale, 0.0f)) || enter_location) {
        navigate_browser(browser, std::filesystem::path(browser.location));
    }
    ImGui::SameLine();
    if (ImGui::Button("Up", ImVec2(64.0f * scale, 0.0f))) {
        const std::filesystem::path parent = browser.current_directory.parent_path();
        if (!parent.empty()) navigate_browser(browser, parent);
    }

    if (!browser.error.empty()) {
        ImGui::TextColored(ImVec4(1.0f, 0.45f, 0.4f, 1.0f), "%s", browser.error.c_str());
    }
    const bool selects_directory = browser.target == PathTarget::InputDirectory ||
                                   browser.target == PathTarget::OutputDirectory;
    if (!selects_directory) {
        ImGui::Checkbox("Show all files", &browser.show_all_files);
        ImGui::SameLine();
        ImGui::TextDisabled(browser.target == PathTarget::FastaFile
                                ? "FASTA files are shown by default"
                                : "raw, .d.zip, and HDF5 files are shown by default");
    }

    std::vector<BrowserEntry> entries;
    std::error_code error;
    for (const auto& entry : std::filesystem::directory_iterator(browser.current_directory, error)) {
        std::error_code entry_error;
        const bool directory = entry.is_directory(entry_error);
        if (entry_error) continue;
        if (!directory && !entry.is_regular_file(entry_error)) continue;
        if (!directory && selects_directory) continue;
        if (!directory && !selects_directory && !browser.show_all_files &&
            !browser_accepts_file(browser.target, entry.path())) continue;
        BrowserEntry item;
        item.path = entry.path();
        item.directory = directory;
        if (!directory) item.size = entry.file_size(entry_error);
        entries.push_back(std::move(item));
    }
    if (error) browser.error = "Unable to read directory: " + error.message();
    std::sort(entries.begin(), entries.end(), [](const BrowserEntry& left,
                                                  const BrowserEntry& right) {
        if (left.directory != right.directory) return left.directory > right.directory;
        return lower_copy(left.path.filename().string()) <
               lower_copy(right.path.filename().string());
    });

    const float footer_height = ImGui::GetFrameHeightWithSpacing() * 2.0f;
    if (ImGui::BeginTable("##browser_entries", 3,
                          ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_RowBg |
                              ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY,
                          ImVec2(0.0f, std::max(150.0f * scale,
                                               ImGui::GetContentRegionAvail().y - footer_height)))) {
        ImGui::TableSetupColumn("Name", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthFixed, 100.0f * scale);
        ImGui::TableSetupColumn("Size", ImGuiTableColumnFlags_WidthFixed, 90.0f * scale);
        ImGui::TableHeadersRow();
        for (const auto& entry : entries) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            const std::string name = (entry.directory ? "[DIR] " : "") +
                                     entry.path.filename().string();
            const std::string id = name + "##" + entry.path.string();
            const bool selected = browser.selected_path == entry.path;
            if (ImGui::Selectable(id.c_str(), selected,
                                  ImGuiSelectableFlags_SpanAllColumns |
                                      ImGuiSelectableFlags_AllowDoubleClick)) {
                browser.selected_path = entry.path;
                if (entry.directory && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    navigate_browser(browser, entry.path);
                } else if (!entry.directory &&
                           ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    assign_path(browser.target, entry.path.string(), form);
                    ImGui::CloseCurrentPopup();
                }
            }
            ImGui::TableSetColumnIndex(1);
            ImGui::TextUnformatted(entry.directory ? "Folder" : "File");
            ImGui::TableSetColumnIndex(2);
            if (!entry.directory) ImGui::TextUnformatted(display_size(entry.size).c_str());
        }
        ImGui::EndTable();
    }

    bool can_select = selects_directory ||
                      (!browser.selected_path.empty() &&
                       std::filesystem::is_regular_file(browser.selected_path, error));
    ImGui::BeginDisabled(!can_select);
    const char* select_label = selects_directory ? "Use this folder" : "Select file";
    if (ImGui::Button(select_label, ImVec2(150.0f * scale, 0.0f))) {
        const std::filesystem::path selection = selects_directory
            ? browser.current_directory
            : browser.selected_path;
        assign_path(browser.target, selection.string(), form);
        ImGui::CloseCurrentPopup();
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button("Cancel", ImVec2(100.0f * scale, 0.0f))) {
        ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
}

void render_mode_selector(WorkflowOptions& form, float scale) {
    ImGui::SeparatorText("Search mode");
    int mode = form.mode == WorkflowMode::RegularFasta ? 0 :
               form.mode == WorkflowMode::FastSip ? 1 : 2;
    ImGui::PushStyleVar(ImGuiStyleVar_SelectableTextAlign, ImVec2(0.5f, 0.5f));
    if (ImGui::BeginTable("##search_modes", 3,
                          ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_BordersInnerV)) {
        const char* labels[] = {
            "Regular FASTA##regular_mode",
            "Fast SIP##fast_mode",
            "SIP FASTA##sip_mode",
        };
        for (int index = 0; index < 3; ++index) {
            ImGui::TableNextColumn();
            if (ImGui::Selectable(labels[index], mode == index, 0,
                                  ImVec2(0.0f, 34.0f * scale))) mode = index;
        }
        ImGui::EndTable();
    }
    ImGui::PopStyleVar();
    form.mode = mode == 0 ? WorkflowMode::RegularFasta :
                mode == 1 ? WorkflowMode::FastSip : WorkflowMode::SipFasta;
}

void render_path_rows(GLFWwindow* window, WorkflowOptions& form,
                      BrowserState& browser, GuiState& state, float scale) {
    ImGui::SeparatorText("Inputs and output");
    if (!ImGui::BeginTable("##workflow_paths", 4,
                           ImGuiTableFlags_SizingStretchProp |
                               ImGuiTableFlags_BordersInnerH |
                               ImGuiTableFlags_PadOuterX)) return;
    const float label_width = std::max(
        170.0f * scale,
        ImGui::CalcTextSize("Raw / HDF5 input").x + 16.0f * scale);
    const auto button_width = [](const char* label) {
        return std::ceil(ImGui::CalcTextSize(label).x +
                         2.0f * ImGui::GetStyle().FramePadding.x + 4.0f);
    };
    const float file_width = std::max(button_width("Files..."),
                                      button_width("Browse..."));
    const float folder_width = std::max(button_width("Folder..."),
                                        button_width("Browse..."));
    ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthFixed, label_width);
    ImGui::TableSetupColumn("Path", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableSetupColumn("File", ImGuiTableColumnFlags_WidthFixed, file_width);
    ImGui::TableSetupColumn("Folder", ImGuiTableColumnFlags_WidthFixed, folder_width);

    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Raw / HDF5 input");
    ImGui::TableSetColumnIndex(1);
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputTextWithHint("##input", "File, directory, or comma-separated files", &form.input);
    ImGui::TableSetColumnIndex(2);
    if (ImGui::Button("Files...##input_files", ImVec2(-FLT_MIN, 0.0f))) {
        request_path(window, PathTarget::InputFiles, form, browser, state);
    }
    ImGui::TableSetColumnIndex(3);
    if (ImGui::Button("Folder...##input_folder", ImVec2(-FLT_MIN, 0.0f))) {
        request_path(window, PathTarget::InputDirectory, form, browser, state);
    }

    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("FASTA database");
    ImGui::TableSetColumnIndex(1);
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputTextWithHint("##fasta", "Target FASTA file", &form.fasta);
    ImGui::TableSetColumnIndex(2);
    if (ImGui::Button("Browse...##fasta_file", ImVec2(-FLT_MIN, 0.0f))) {
        request_path(window, PathTarget::FastaFile, form, browser, state);
    }

    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Output directory");
    ImGui::TableSetColumnIndex(1);
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputTextWithHint("##output", "Workflow output directory", &form.output);
    ImGui::TableSetColumnIndex(3);
    if (ImGui::Button("Browse...##output_folder", ImVec2(-FLT_MIN, 0.0f))) {
        request_path(window, PathTarget::OutputDirectory, form, browser, state);
    }
    ImGui::EndTable();
}

void render_sip_options(WorkflowOptions& form, int& isotope_index,
                        std::string& sip_range, std::string& precision,
                        const char* const* isotopes, int isotope_count,
                        float scale) {
    if (form.mode == WorkflowMode::RegularFasta) {
        form.element = "R";
        form.sip_range.reset();
        form.precision.reset();
        return;
    }
    ImGui::SeparatorText("SIP labeling");
    if (ImGui::BeginTable("##sip_options", 3, ImGuiTableFlags_SizingStretchSame)) {
        ImGui::TableNextColumn();
        ImGui::TextUnformatted("Isotope");
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::Combo("##sip_isotope", &isotope_index, isotopes, isotope_count);
        ImGui::TableNextColumn();
        ImGui::TextUnformatted("Label range (%)");
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::InputText("##sip_range", &sip_range);
        ImGui::TableNextColumn();
        ImGui::TextUnformatted("Precision (%)");
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::InputText("##sip_precision", &precision);
        ImGui::EndTable();
    }
    (void)scale;
    form.element = isotopes[isotope_index];
    form.sip_range = sip_range.empty() ? std::nullopt : std::optional<std::string>(sip_range);
    form.precision = precision.empty() ? std::nullopt : std::optional<std::string>(precision);
}

void render_advanced_options(GLFWwindow* window, WorkflowOptions& form,
                             PtmSelection& ptm_selection, int& max_ptm_count,
                             BrowserState& browser, GuiState& state,
                             float scale) {
    if (!ImGui::CollapsingHeader("Advanced options")) return;
    const int columns = ImGui::GetContentRegionAvail().x >= 760.0f * scale ? 2 : 1;
    if (ImGui::BeginTable("##advanced_values", columns,
                          ImGuiTableFlags_SizingStretchSame |
                              (columns == 2 ? ImGuiTableFlags_BordersInnerV : 0))) {
        auto input_double = [&](const char* label, const char* id, double* value,
                                double step, double fast_step,
                                const char* format) {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(label);
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::InputDouble(id, value, step, fast_step, format);
        };
        auto input_int = [&](const char* label, const char* id, int* value) {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(label);
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::InputInt(id, value);
        };
        input_double("MS1 tolerance (Da)", "##tolerance_ms1", &form.tolerance_ms1,
                     0.001, 0.01, "%.6f");
        input_double("MS2 tolerance (Da)", "##tolerance_ms2", &form.tolerance_ms2,
                     0.001, 0.01, "%.6f");
        input_int("Physical-core budget", "##threads", &form.threads);
        input_int("Raxport precursors", "##n_precursor", &form.n_precursor);
        input_int("Top PSMs per scan", "##top_psms", &form.top_psms_per_scan);
        input_int("Product top isotopes", "##product_top", &form.product_top_isotopes);
        input_int("Aerith sample parallelism", "##aerith_parallel",
                  &form.aerith_sample_parallelism);
        input_double("RT tolerance (min)", "##rt_tolerance", &form.rt_tolerance,
                     0.5, 1.0, "%.2f");
        input_int("SFI envelope top N", "##sfi_top", &form.sfi_envelope_top_n);
        input_int("MVH cascade top N", "##mvh_top", &form.mvh_cascade_top_n);
        if (form.mode != WorkflowMode::RegularFasta) {
            input_double("Label threshold (%)", "##label_threshold",
                         &form.label_threshold, 0.1, 1.0, "%.2f");
            ImGui::TableNextColumn();
            ImGui::Dummy(ImVec2(0.0f, ImGui::GetTextLineHeight()));
            ImGui::Checkbox("Ignore SIP abundance", &form.ignore_pct);
        }
        ImGui::EndTable();
    }

    if (form.mode != WorkflowMode::RegularFasta) {
        ImGui::SeparatorText("Negative controls");
        if (ImGui::BeginTable("##negative_control_picker", 2,
                              ImGuiTableFlags_SizingStretchProp)) {
            ImGui::TableSetupColumn("Samples", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableSetupColumn("Picker", ImGuiTableColumnFlags_WidthFixed,
                                    150.0f * scale);
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::InputTextWithHint("##negative", "Selected sample basenames",
                                     &form.negative_control);
            ImGui::TableNextColumn();
            if (ImGui::Button("Select files...##negative_files",
                              ImVec2(-FLT_MIN, 0.0f))) {
                request_path(window, PathTarget::NegativeControlFiles, form,
                             browser, state);
            }
            ImGui::EndTable();
        }
        ImGui::TextDisabled("Select one or more raw/HDF5 files; their sample basenames are used.");
    }

    if (form.mode == WorkflowMode::RegularFasta) {
        ImGui::SeparatorText("Regular-search PTMs");
        if (ImGui::Checkbox("Fixed carbamidomethylation (C)",
                            &ptm_selection.fixed_carbamidomethyl) &&
            ptm_selection.fixed_carbamidomethyl) {
            ptm_selection.variable[11] = false;
        }
        ImGui::TextDisabled("Select any compatible variable PTMs:");
        const float available = ImGui::GetContentRegionAvail().x;
        const int ptm_columns = available >= 900.0f * scale ? 3 :
                                available >= 560.0f * scale ? 2 : 1;
        if (ImGui::BeginTable("##regular_ptm_choices", ptm_columns,
                              ImGuiTableFlags_SizingStretchSame)) {
            for (std::size_t index = 0; index < kVariablePtmChoices.size(); ++index) {
                ImGui::TableNextColumn();
                const bool incompatible_iaa = index == 11 &&
                    ptm_selection.fixed_carbamidomethyl;
                ImGui::BeginDisabled(incompatible_iaa);
                ImGui::Checkbox(kVariablePtmChoices[index].label,
                                &ptm_selection.variable[index]);
                ImGui::EndDisabled();
            }
            ImGui::EndTable();
        }
        ImGui::Text("Maximum variable PTMs per peptide (default: %d)",
                    kRegularDefaultMaxPtmCount);
        ImGui::SetNextItemWidth(180.0f * scale);
        ImGui::InputInt("##max_ptm_count", &max_ptm_count);
        max_ptm_count = std::max(0, max_ptm_count);
    }

    ImGui::Checkbox("Dry run", &form.dry_run);
}

void sync_ptm_options(WorkflowOptions& form, const PtmSelection& selection,
                      int max_ptm_count) {
    form.ptms.clear();
    form.fixed_ptms.clear();
    form.max_ptm_count.reset();
    if (form.mode != WorkflowMode::RegularFasta) return;

    for (std::size_t index = 0; index < kVariablePtmChoices.size(); ++index) {
        if (selection.variable[index]) {
            form.ptms.emplace_back(kVariablePtmChoices[index].selector);
        }
    }
    if (form.ptms.empty()) form.ptms.emplace_back("none");
    form.fixed_ptms.emplace_back(
        selection.fixed_carbamidomethyl ? "carbamidomethyl" : "none");
    form.max_ptm_count = max_ptm_count;
}

void reset_form(WorkflowOptions& form, std::string& sip_range,
                std::string& precision, PtmSelection& ptm_selection,
                int& max_ptm_count, int& isotope_index) {
    form = WorkflowOptions{};
    form.mode = WorkflowMode::RegularFasta;
    form.element = "C13";
    sip_range = "0-100";
    precision = "1";
    ptm_selection = PtmSelection{};
    ptm_selection.variable[0] = true;
    max_ptm_count = kRegularDefaultMaxPtmCount;
    isotope_index = 0;
}

void launch_workflow(GuiState& state, WorkflowOptions options,
                     const std::filesystem::path& executable_path) {
    std::vector<std::string> warnings;
    try {
        if (options.mode == WorkflowMode::RegularFasta) {
            options.negative_control.clear();
            options.ignore_pct = false;
        }
        validate_options(options, &warnings);
        state.cancelled.store(false);
        state.running.store(true);
        set_status(state, "Running " + mode_name(options.mode));
        state.worker = std::thread([&state, options = std::move(options), warnings,
                                    executable_path] {
            try {
                std::error_code error;
                const bool existed = std::filesystem::exists(options.output, error);
                std::filesystem::create_directories(options.output);
                const std::filesystem::path log_path =
                    std::filesystem::path(options.output) / kWorkflowLogFilename;
                Logger logger(log_path, gui_log_sink(state));
                try {
                    logger.info("Workflow log: " + log_path.string());
                    if (existed) logger.warning(options.output + " exists and will be overwritten");
                    for (const auto& warning : warnings) logger.warning(warning);
                    Workflow workflow(options, locate_tools(executable_path), logger,
                                      state.cancelled);
                    workflow.run();
                    set_status(state, "Completed successfully");
                } catch (const std::exception& error) {
                    if (state.cancelled.load()) {
                        logger.warning(std::string("Workflow cancelled: ") + error.what());
                        set_status(state, "Cancelled");
                    } else {
                        logger.error(std::string("Workflow failed: ") + error.what());
                        set_status(state, "Failed");
                    }
                }
            } catch (const std::exception& error) {
                append_log(state, std::string("ERROR: ") + error.what());
                set_status(state, "Failed");
            }
            state.running.store(false);
            glfwPostEmptyEvent();
        });
    } catch (const std::exception& error) {
        state.running.store(false);
        append_log(state, std::string("Validation error: ") + error.what());
        set_status(state, "Validation failed");
    }
}

bool render_command_and_actions(WorkflowOptions& form, GuiState& state,
                                const std::filesystem::path& executable_path, float scale) {
    ImGui::SeparatorText("Command preview");
    const float preview_height = ImGui::GetTextLineHeightWithSpacing() * 3.7f;
    if (ImGui::BeginChild("##command_preview", ImVec2(0.0f, preview_height),
                          ImGuiChildFlags_Borders | ImGuiChildFlags_AlwaysUseWindowPadding)) {
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextUnformatted(preview(form).c_str());
        ImGui::PopTextWrapPos();
    }
    ImGui::EndChild();

    const bool running = state.running.load();
    bool reset_requested = false;
    ImGui::BeginDisabled(running);
    if (ImGui::Button("Run workflow", ImVec2(160.0f * scale, 0.0f))) {
        launch_workflow(state, form, executable_path);
    }
    ImGui::SameLine();
    if (ImGui::Button("Reset defaults", ImVec2(170.0f * scale, 0.0f))) {
        reset_requested = true;
    }
    ImGui::EndDisabled();
    if (running) {
        ImGui::SameLine();
        if (ImGui::Button("Stop", ImVec2(100.0f * scale, 0.0f))) {
            state.cancelled.store(true);
            set_status(state, "Stopping active processes...");
        }
    }
    ImGui::SameLine();
    const std::string status = current_status(state);
    ImVec4 color = ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled);
    if (status.find("Completed") != std::string::npos) color = ImVec4(0.35f, 0.90f, 0.55f, 1.0f);
    else if (status.find("Failed") != std::string::npos ||
             status.find("error") != std::string::npos) color = ImVec4(1.0f, 0.45f, 0.4f, 1.0f);
    else if (running) color = ImVec4(0.4f, 0.72f, 1.0f, 1.0f);
    ImGui::TextColored(color, "%s", status.c_str());
    return reset_requested;
}

void render_log(GuiState& state, float scale) {
    std::size_t shown_lines = 0;
    std::size_t hidden_lines = 0;
    {
        std::lock_guard<std::mutex> lock(state.mutex);
        shown_lines = state.logs.size();
        hidden_lines = state.hidden_log_lines;
    }

    ImGui::Spacing();
    if (ImGui::BeginTable("##log_header", 2, ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Title", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Action", ImGuiTableColumnFlags_WidthFixed, 92.0f * scale);
        ImGui::TableNextColumn();
        ImGui::AlignTextToFramePadding();
        if (hidden_lines == 0) {
            ImGui::Text("Workflow log  (%zu lines)", shown_lines);
        } else {
            ImGui::Text("Workflow log  (%zu shown, %zu older hidden)",
                        shown_lines, hidden_lines);
        }
        ImGui::TableNextColumn();
        if (ImGui::Button("Clear##clear_log", ImVec2(-FLT_MIN, 0.0f))) {
            std::lock_guard<std::mutex> lock(state.mutex);
            state.logs.clear();
            state.hidden_log_lines = 0;
            state.scroll_to_bottom = false;
        }
        ImGui::EndTable();
    }
    ImGui::Separator();
    const float height = std::max(180.0f * scale, ImGui::GetContentRegionAvail().y);
    if (ImGui::BeginChild("##workflow_log", ImVec2(0.0f, height),
                          ImGuiChildFlags_Borders | ImGuiChildFlags_AlwaysUseWindowPadding,
                          ImGuiWindowFlags_HorizontalScrollbar)) {
        std::lock_guard<std::mutex> lock(state.mutex);
        if (state.logs.empty()) ImGui::TextDisabled("Workflow output will appear here.");
        else {
            ImGuiListClipper clipper;
            clipper.Begin(static_cast<int>(state.logs.size()));
            while (clipper.Step()) {
                for (int index = clipper.DisplayStart;
                     index < clipper.DisplayEnd; ++index) {
                    const std::string& line =
                        state.logs[static_cast<std::size_t>(index)];
                    ImVec4 color = ImGui::GetStyleColorVec4(ImGuiCol_Text);
                    if (line.find(" | ERROR | ") != std::string::npos ||
                        line.rfind("ERROR:", 0) == 0) {
                        color = ImVec4(1.0f, 0.43f, 0.38f, 1.0f);
                    } else if (line.find(" | WARN  | ") != std::string::npos) {
                        color = ImVec4(1.0f, 0.76f, 0.32f, 1.0f);
                    } else if (line.find(" | DEBUG | ") != std::string::npos) {
                        color = ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled);
                    }
                    ImGui::PushStyleColor(ImGuiCol_Text, color);
                    ImGui::TextUnformatted(line.c_str());
                    ImGui::PopStyleColor();
                }
            }
        }
        if (state.scroll_to_bottom) {
            ImGui::SetScrollHereY(1.0f);
            state.scroll_to_bottom = false;
        }
    }
    ImGui::EndChild();
}

MonitorWorkArea monitor_work_area(GLFWmonitor* monitor) {
    MonitorWorkArea area;
    if (monitor != nullptr) {
        glfwGetMonitorWorkarea(monitor, &area.x, &area.y, &area.width,
                               &area.height);
        if (area.width <= 0 || area.height <= 0) area = MonitorWorkArea{};
    }
    return area;
}

GLFWmonitor* monitor_for_window(GLFWwindow* window) {
    if (GLFWmonitor* monitor = glfwGetWindowMonitor(window); monitor != nullptr) {
        return monitor;
    }
#if defined(__linux__) && defined(GLFW_PLATFORM_WAYLAND)
    if (glfwGetPlatform() == GLFW_PLATFORM_WAYLAND) {
        return glfwGetPrimaryMonitor();
    }
#endif
    int window_x = 0, window_y = 0, window_width = 0, window_height = 0;
    glfwGetWindowPos(window, &window_x, &window_y);
    glfwGetWindowSize(window, &window_width, &window_height);

    int monitor_count = 0;
    GLFWmonitor** monitors = glfwGetMonitors(&monitor_count);
    GLFWmonitor* best_monitor = nullptr;
    std::int64_t best_overlap = -1;
    for (int index = 0; index < monitor_count; ++index) {
        const MonitorWorkArea area = monitor_work_area(monitors[index]);
        const int overlap_width = std::max(
            0, std::min(window_x + window_width, area.x + area.width) -
                   std::max(window_x, area.x));
        const int overlap_height = std::max(
            0, std::min(window_y + window_height, area.y + area.height) -
                   std::max(window_y, area.y));
        const std::int64_t overlap =
            static_cast<std::int64_t>(overlap_width) * overlap_height;
        if (overlap > best_overlap) {
            best_overlap = overlap;
            best_monitor = monitors[index];
        }
    }
    return best_monitor != nullptr ? best_monitor : glfwGetPrimaryMonitor();
}

float fitted_content_scale(float content_scale, const MonitorWorkArea& area) {
    const float normalized = normalized_content_scale(content_scale, content_scale);
    const float width_fit =
        static_cast<float>(std::max(640, area.width - 48)) / 1024.0f;
    const float height_fit =
        static_cast<float>(std::max(480, area.height - 48)) / 768.0f;
    return std::clamp(std::min({normalized, width_fit, height_fit}), 1.0f,
                      normalized);
}

void configure_default_font() {
    ImGuiIO& io = ImGui::GetIO();
    io.FontDefault = io.Fonts->AddFontDefaultVector();
}

void configure_style(float scale) {
    ImGuiStyle& style = ImGui::GetStyle();
    style = ImGuiStyle{};
    ImGui::StyleColorsDark();
    style.ScaleAllSizes(scale);
    style.FontSizeBase = kBaseFontSize;
    style.FontScaleDpi = scale;
    style.WindowRounding = 0.0f;
    style.ChildRounding = 5.0f * scale;
    style.FrameRounding = 4.0f * scale;
    style.PopupRounding = 6.0f * scale;
    style.ScrollbarRounding = 5.0f * scale;
    style.GrabRounding = 4.0f * scale;
    style.WindowBorderSize = 0.0f;
    style.FrameBorderSize = 1.0f;
    style.FramePadding = ImVec2(8.0f * scale, 5.0f * scale);
    style.ItemSpacing = ImVec2(8.0f * scale, 7.0f * scale);
    style.WindowPadding = ImVec2(14.0f * scale, 12.0f * scale);
    style.Colors[ImGuiCol_Text] = ImVec4(0.91f, 0.96f, 0.97f, 1.0f);
    style.Colors[ImGuiCol_TextDisabled] = ImVec4(0.54f, 0.63f, 0.65f, 1.0f);
    style.Colors[ImGuiCol_WindowBg] = ImVec4(0.035f, 0.052f, 0.058f, 1.0f);
    style.Colors[ImGuiCol_ChildBg] = ImVec4(0.043f, 0.068f, 0.075f, 1.0f);
    style.Colors[ImGuiCol_PopupBg] = ImVec4(0.050f, 0.078f, 0.086f, 0.98f);
    style.Colors[ImGuiCol_TitleBg] = ImVec4(0.055f, 0.090f, 0.100f, 1.0f);
    style.Colors[ImGuiCol_TitleBgActive] = ImVec4(0.145f, 0.315f, 0.340f, 1.0f);
    style.Colors[ImGuiCol_TitleBgCollapsed] = ImVec4(0.055f, 0.090f, 0.100f, 1.0f);
    style.Colors[ImGuiCol_FrameBg] = ImVec4(0.075f, 0.130f, 0.142f, 1.0f);
    style.Colors[ImGuiCol_FrameBgHovered] = ImVec4(0.115f, 0.225f, 0.240f, 1.0f);
    style.Colors[ImGuiCol_FrameBgActive] = ImVec4(0.145f, 0.285f, 0.305f, 1.0f);
    style.Colors[ImGuiCol_Button] = ImVec4(0.145f, 0.315f, 0.340f, 1.0f);
    style.Colors[ImGuiCol_ButtonHovered] = ImVec4(0.205f, 0.455f, 0.485f, 1.0f);
    style.Colors[ImGuiCol_ButtonActive] = ImVec4(0.245f, 0.545f, 0.575f, 1.0f);
    style.Colors[ImGuiCol_Header] = ImVec4(0.165f, 0.365f, 0.390f, 1.0f);
    style.Colors[ImGuiCol_HeaderHovered] = ImVec4(0.225f, 0.500f, 0.525f, 1.0f);
    style.Colors[ImGuiCol_HeaderActive] = ImVec4(0.265f, 0.575f, 0.600f, 1.0f);
    style.Colors[ImGuiCol_CheckMark] = ImVec4(0.62f, 0.90f, 0.92f, 1.0f);
    style.Colors[ImGuiCol_SliderGrab] = ImVec4(0.48f, 0.73f, 0.76f, 1.0f);
    style.Colors[ImGuiCol_SliderGrabActive] = ImVec4(0.65f, 0.91f, 0.92f, 1.0f);
    style.Colors[ImGuiCol_Border] = ImVec4(0.40f, 0.55f, 0.58f, 0.75f);
    style.Colors[ImGuiCol_Separator] = ImVec4(0.34f, 0.48f, 0.51f, 0.80f);
    style.Colors[ImGuiCol_ResizeGrip] = ImVec4(0.42f, 0.67f, 0.70f, 0.55f);
    style.Colors[ImGuiCol_ResizeGripHovered] = ImVec4(0.60f, 0.88f, 0.90f, 0.85f);
}

void update_window_size_limits(GLFWwindow* window, float scale,
                               int& minimum_width, int& minimum_height) {
    const MonitorWorkArea area = monitor_work_area(monitor_for_window(window));
    minimum_width = std::min(
        area.width, static_cast<int>(std::lround(kBaseMinimumWidth * scale)));
    minimum_height = std::min(
        area.height, static_cast<int>(std::lround(kBaseMinimumHeight * scale)));
    glfwSetWindowSizeLimits(window, minimum_width, minimum_height,
                            GLFW_DONT_CARE, GLFW_DONT_CARE);
}

bool update_runtime_dpi_scale(GLFWwindow* window, WindowState& window_state,
                              float& current_scale,
                              unsigned int& title_icon_texture,
                              int& minimum_width, int& minimum_height) {
    float scale_x = 1.0f;
    float scale_y = 1.0f;
    glfwGetWindowContentScale(window, &scale_x, &scale_y);
    const float observed_scale = normalized_content_scale(scale_x, scale_y);
    if (std::abs(observed_scale - window_state.content_scale) >
        kDpiScaleEpsilon) {
        window_state.content_scale = observed_scale;
        window_state.dpi_change_pending = true;
    }

    const MonitorWorkArea area = monitor_work_area(monitor_for_window(window));
    const float target_scale = fitted_content_scale(
        window_state.content_scale, area);
    if (!window_state.dpi_change_pending &&
        std::abs(target_scale - current_scale) <= kDpiScaleEpsilon) {
        return false;
    }
    window_state.dpi_change_pending = false;
    if (std::abs(target_scale - current_scale) <= kDpiScaleEpsilon) {
        return false;
    }

    current_scale = target_scale;
    configure_style(current_scale);
    update_window_size_limits(window, current_scale, minimum_width,
                              minimum_height);
    if (title_icon_texture != 0) {
        glDeleteTextures(1, &title_icon_texture);
    }
    title_icon_texture = create_title_icon_texture(
        std::max(32, static_cast<int>(std::lround(64.0f * current_scale))));
    install_native_window_icon(window);
    return true;
}

int run_gui(const std::filesystem::path& executable_path) {
    bool use_custom_frame = true;
    if (!initialize_glfw(use_custom_frame)) {
        throw std::runtime_error("Unable to initialize GLFW");
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_CONTEXT_CREATION_API, GLFW_NATIVE_CONTEXT_API);
    glfwWindowHint(GLFW_DECORATED,
                   use_custom_frame ? GLFW_FALSE : GLFW_TRUE);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);

    GLFWmonitor* monitor = glfwGetPrimaryMonitor();
    const MonitorWorkArea initial_work_area = monitor_work_area(monitor);
    float monitor_scale_x = 1.0f;
    float monitor_scale_y = 1.0f;
    if (monitor != nullptr) {
        glfwGetMonitorContentScale(monitor, &monitor_scale_x, &monitor_scale_y);
    }
    const float system_scale = normalized_content_scale(
        monitor_scale_x, monitor_scale_y);
    const float initial_scale = fitted_content_scale(system_scale,
                                                       initial_work_area);
    const int margin = static_cast<int>(std::lround(24.0f * initial_scale));
    const int desired_width = static_cast<int>(std::lround(1024.0f * initial_scale));
    const int desired_height = static_cast<int>(std::lround(768.0f * initial_scale));
    const int window_width = std::max(
        640, std::min(desired_width, initial_work_area.width - margin));
    const int window_height = std::max(
        480, std::min(desired_height, initial_work_area.height - margin));
    GLFWwindow* window = glfwCreateWindow(window_width, window_height,
                                          "Sipros Workflow", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        throw std::runtime_error("Unable to create the Sipros Workflow window");
    }
    install_native_window_icon(window);
    if (use_custom_frame) {
        glfwSetWindowPos(
            window,
            initial_work_area.x +
                std::max(0, (initial_work_area.width - window_width) / 2),
            initial_work_area.y +
                std::max(0, (initial_work_area.height - window_height) / 2));
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    WindowState window_state;
    float current_scale_x = 1.0f;
    float current_scale_y = 1.0f;
    glfwGetWindowContentScale(window, &current_scale_x, &current_scale_y);
    window_state.content_scale = normalized_content_scale(
        current_scale_x, current_scale_y);
    float current_scale = fitted_content_scale(
        window_state.content_scale,
        monitor_work_area(monitor_for_window(window)));
    int minimum_width = 0;
    int minimum_height = 0;
    update_window_size_limits(window, current_scale, minimum_width,
                              minimum_height);
    unsigned int title_icon_texture = create_title_icon_texture(
        std::max(32, static_cast<int>(std::lround(64.0f * current_scale))));

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.IniFilename = nullptr;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    configure_default_font();
    configure_style(current_scale);
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    WorkflowOptions form;
    std::string sip_range;
    std::string precision;
    PtmSelection ptm_selection;
    int max_ptm_count = kRegularDefaultMaxPtmCount;
    int isotope_index = 0;
    const char* isotopes[] = {"C13", "H2", "N15", "O18", "S34"};
    GuiState state;
    BrowserState browser;
    FramelessState frameless;
    reset_form(form, sip_range, precision, ptm_selection,
               max_ptm_count, isotope_index);
    glfwSetWindowUserPointer(window, &window_state);
    glfwSetWindowCloseCallback(window, handle_window_close);
    window_state.previous_content_scale_callback =
        glfwSetWindowContentScaleCallback(window,
                                          handle_window_content_scale);

    bool first_frame = true;
    while (!glfwWindowShouldClose(window)) {
        if (first_frame) {
            glfwPollEvents();
            first_frame = false;
        } else {
            const bool iconified =
                glfwGetWindowAttrib(window, GLFW_ICONIFIED) == GLFW_TRUE;
            if (state.running.load() && !iconified) {
                glfwWaitEventsTimeout(kActiveEventWaitSeconds);
            } else {
                glfwWaitEvents();
            }
        }
        if (glfwGetWindowAttrib(window, GLFW_ICONIFIED) == GLFW_TRUE) continue;
        if (update_runtime_dpi_scale(
                window, window_state, current_scale, title_icon_texture,
                minimum_width, minimum_height)) {
            frameless = FramelessState{};
        }
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        if (!state.running.load() && state.worker.joinable()) state.worker.join();
        const float title_height = use_custom_frame
            ? kBaseTitleHeight * current_scale
            : 0.0f;
        const float title_button_width = kBaseTitleButtonWidth * current_scale;
        const float title_controls_width = 3.0f * title_button_width;
        int resize_cursor = ResizeNone;
        if (use_custom_frame) {
            resize_cursor = update_frameless_interaction(
                window, frameless, title_height, title_controls_width,
                kBaseResizeBorder * current_scale, minimum_width,
                minimum_height);
        }

        ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f));
        ImGui::SetNextWindowSize(io.DisplaySize);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
        ImGui::Begin("##sipros_root", nullptr,
                     ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
                         ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings |
                         ImGuiWindowFlags_NoBringToFrontOnFocus);
        ImGui::PopStyleVar();
        if (use_custom_frame) {
            render_title_bar(window, current_scale, title_height,
                             title_button_width, title_icon_texture,
                             window_state.exit);
        }

        ImGui::SetCursorPos(ImVec2(0.0f, title_height));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding,
                            ImVec2(14.0f * current_scale, 12.0f * current_scale));
        const bool content_visible = ImGui::BeginChild(
            "##workflow_content", ImVec2(0.0f, 0.0f),
            ImGuiChildFlags_AlwaysUseWindowPadding, ImGuiWindowFlags_None);
        ImGui::PopStyleVar();
        if (content_visible) {
            const bool running = state.running.load();
            ImGui::BeginDisabled(running);
            render_mode_selector(form, current_scale);
            render_path_rows(window, form, browser, state, current_scale);
            render_sip_options(form, isotope_index, sip_range, precision,
                               isotopes, IM_ARRAYSIZE(isotopes), current_scale);
            render_advanced_options(window, form, ptm_selection, max_ptm_count,
                                    browser, state, current_scale);
            ImGui::EndDisabled();
            sync_ptm_options(form, ptm_selection, max_ptm_count);
            if (render_command_and_actions(form, state, executable_path,
                                           current_scale)) {
                reset_form(form, sip_range, precision, ptm_selection,
                           max_ptm_count, isotope_index);
                set_status(state, "Ready");
            }
            render_log(state, current_scale);
        }
        ImGui::EndChild();

        if (use_custom_frame) {
            const ImVec2 window_size = io.DisplaySize;
            ImGui::GetWindowDrawList()->AddTriangleFilled(
                ImVec2(window_size.x, window_size.y),
                ImVec2(window_size.x - 13.0f * current_scale, window_size.y),
                ImVec2(window_size.x, window_size.y - 13.0f * current_scale),
                ImGui::GetColorU32(ImGuiCol_ResizeGrip));
        }
        ImGui::End();

        render_browser(browser, form, current_scale);
        render_exit_confirmation(window, window_state.exit, state,
                                 current_scale);
        if (resize_cursor != ResizeNone) set_resize_cursor(resize_cursor);

        ImGui::Render();
        int framebuffer_width = 0, framebuffer_height = 0;
        glfwGetFramebufferSize(window, &framebuffer_width, &framebuffer_height);
        glViewport(0, 0, framebuffer_width, framebuffer_height);
        glClearColor(0.035f, 0.052f, 0.058f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    state.cancelled.store(true);
    if (state.worker.joinable()) state.worker.join();
    glfwSetWindowContentScaleCallback(
        window, window_state.previous_content_scale_callback);
    glfwSetWindowCloseCallback(window, nullptr);
    glfwSetWindowUserPointer(window, nullptr);
    if (title_icon_texture != 0) glDeleteTextures(1, &title_icon_texture);
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

} // namespace siproswf
