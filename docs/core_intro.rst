Getting started
====================

PSDR-CUDA is a physics-based differenetiable renderer written in C++17/CUDA and based on the theoretical framework introduced in the paper titled "`Path-Space Differentiable Rendering <https://shuangz.com/projects/psdr-sg20/>`_".

PSDR-CUDA is light-weight (compared to other modern renderers like `Mitsuba 2 <https://www.mitsuba-renderer.org/>`_) yet powerful, offering the following features:

\1. **Fast GPU-based ray tracing:** PSDR-CUDA uses `OptiX 7.1 <https://developer.nvidia.com/optix/>`_ that leverages RTX ray tracing on modern Nvidia graphics hardware.

\2. **Vectorized computation:** PSDR-CUDA uses `Enoki <https://github.com/mitsuba-renderer/enoki/>`_, the numerical foundation of `Mitsuba 2 <https://www.mitsuba-renderer.org/>`_, to perform vectorized computations (including reverse-mode automatic differentiation) on the GPU using `CUDA 11 <https://developer.nvidia.com/cuda-toolkit>`_.

\3. **Fast and unbiased gradients:** PSDR-CUDA implements state-of-the-art algorithms that produce unbiased gradient estimates with respect to *arbitrary* scene parameters (e.g., spatially varying reflectance and object geometry).

\4. **Python bindings:** PSDR-CUDA provides fine-grained Python bindings using `pybind11 <https://github.com/pybind/pybind11/>`_ and can be easily integrated into Python-based interference and learning pipelines.


About
--------------------

This project was created by `Shuang Zhao <https://shuangz.com/>`_.
Significant features and/or improvements to the code were contributed by Kai Yan.
