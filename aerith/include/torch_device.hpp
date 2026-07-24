#pragma once

#ifdef AERITH_WITH_TORCH

#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <torch/cuda.h>
#include <torch/script.h>

namespace aerith {

inline void compact_diann_gru_weights(torch::jit::script::Module& model) {
    if (!model.hasattr("rnn")) {
        return;
    }
    auto rnn = model.attr("rnn").toModule();
    if (!rnn.hasattr("mode") || rnn.attr("mode").toStringRef() != "GRU" ||
        !rnn.hasattr("_flat_weights") ||
        !at::_use_cudnn_rnn_flatten_weight()) {
        return;
    }

    std::vector<torch::Tensor> weights;
    const auto flat_weights = rnn.attr("_flat_weights").toTensorList();
    weights.reserve(flat_weights.size());
    for (const auto& weight : flat_weights) {
        weights.push_back(weight);
    }
    if (weights.empty()) {
        return;
    }

    // PyTorch's RNNBase::flatten_parameters operation, specialized for the
    // bidirectional GRUs stored in the DIA-NN 2.6.1 TorchScript checkpoints.
    // This packs the weights once instead of having cuDNN repack every batch.
    torch::NoGradGuard no_grad;
    constexpr std::int64_t kGruMode = 3;
    constexpr std::int64_t kWeightsPerLayerAndDirection = 4;
    at::_cudnn_rnn_flatten_weight(
        weights, kWeightsPerLayerAndDirection,
        rnn.attr("input_size").toInt(), kGruMode,
        rnn.attr("hidden_size").toInt(), rnn.attr("proj_size").toInt(),
        rnn.attr("num_layers").toInt(), rnn.attr("batch_first").toBool(),
        rnn.attr("bidirectional").toBool());
}

inline torch::jit::script::Module load_torch_model_on_device(
    const std::filesystem::path& model_path, const torch::Device& device) {
    auto model = torch::jit::load(model_path.string());
    model.eval();
    // DIA-NN's graphs allocate their recurrent hidden state from this string
    // attribute, so moving parameters alone is insufficient for CUDA.
    if (model.hasattr("device")) {
        model.setattr("device", c10::IValue(device.str()));
    }
    model.to(device);
    if (device.is_cuda()) {
        compact_diann_gru_weights(model);
    }
    return model;
}

template <typename Function>
auto run_torch_prefer_cuda(const char* operation, std::string& selected_device,
                           Function&& function)
    -> decltype(function(torch::Device(torch::kCPU))) {
    bool cuda_available = false;
    try {
        cuda_available = torch::cuda::is_available();
    } catch (const std::exception& error) {
        std::cerr << operation << ": CUDA probe failed (" << error.what()
                  << "); using CPU\n";
    }
    if (cuda_available) {
        try {
            std::cerr << operation << ": using CUDA\n";
            auto result = function(torch::Device(torch::kCUDA));
            selected_device = "GPU";
            return result;
        } catch (const std::exception& error) {
            std::cerr << operation << ": CUDA failed (" << error.what()
                      << "); retrying on CPU\n";
        }
    } else {
        std::cerr << operation << ": CUDA unavailable; using CPU\n";
    }
    auto result = function(torch::Device(torch::kCPU));
    selected_device = "CPU";
    return result;
}

inline torch::Tensor move_torch_input(
    const torch::Tensor& cpu_input, const torch::Device& device) {
    return device.is_cpu() ? cpu_input : cpu_input.to(device);
}

} // namespace aerith

#endif
