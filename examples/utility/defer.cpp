#include <vector>

#include "libfly/utility/core.hpp"

void add_count(std::vector<int>& counts) {
  //
  bool commit = false;

  counts.push_back(42);  // (1) direct action.

  // Lambda executed when the enclosing scope exits (function returns).
  fly::Defer _ = [&]() noexcept {
    if (!commit) {
      counts.pop_back();  // (2) rollback action.
    }
  };

  // ...                 // (3) other operations that may throw.

  commit = true;  // ...    (4) disable rollback actions if no throw.
}