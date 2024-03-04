#pragma once

#include <exception>

namespace fly {

  struct not_at_min : std::exception {
    const char* what() const noexcept override { return "Negative eigen-mode"; }
  };

}  // namespace fly