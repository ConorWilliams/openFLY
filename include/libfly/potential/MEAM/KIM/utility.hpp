//
// utility.hpp
//
// LGPL Version 2.1 HEADER START
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
//
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301  USA
//
// LGPL Version 2.1 HEADER END
//

//
// Copyright (c) 2020--2021, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Yaser Afshar
//

#ifndef UTILITY_HPP
#define UTILITY_HPP

#pragma GCC system_header

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

/*!
 * \file utility.hpp
 *
 * \brief The Utility class contains helper functions to read and
 *        parse a file
 *
 */

namespace openKIM {

  /*! \class Utility
   * \brief Utility helper class
   *
   */
  class Utility {
  public:
    /*!
     * \brief Get the first Line.
     *
     * This function does not skip the commented or blank first line
     *
     * \param file_pointer Pointer to the current open file
     * \param next_line_pointer Pointer to the line
     * \param max_line_size Maximum size of the line
     * \return int Flag to indicate if we reached the end of the file
     */
    inline int GetFirstLine(std::FILE *const file_pointer, char *next_line_pointer, int const max_line_size);

    /*!
     * \brief Get the Next Line.
     *
     * It skips empty lines and any commented lines which are lines having a
     * \c # as the first character.
     *
     * \param file_pointer Pointer to the current open file
     * \param next_line_pointer Pointer to the line
     * \param max_line_size Maximum size of the line
     * \return int Flag to indicate if we reached the end of the file
     */
    inline int GetNextLine(std::FILE *const file_pointer, char *next_line_pointer, int const max_line_size);

    /*!
     * \brief Get the number of words in a line.
     *
     * A word is seperated with a white space, or any of the white characters from
     * the reset of the words in a line.
     *
     * \param line_pointer Pointer to the line
     * \return int
     */
    inline int GetNumberOfWordsInLine(char const *line_pointer);

    /*!
     * \brief Get the words in a line.
     *
     * \param buf Pointer to the line
     * \return std::vector<std::string>
     */
    inline std::vector<std::string> GetWordsInLine(char const *line_pointer);
  };

  inline int Utility::GetFirstLine(std::FILE *const file_pointer,
                                   char *next_line_pointer,
                                   int const max_line_size) {
    if (!std::fgets(next_line_pointer, max_line_size, file_pointer)) {
      return 1;
    }
    return 0;
  }

  inline int Utility::GetNextLine(std::FILE *const file_pointer,
                                  char *next_line_pointer,
                                  int const max_line_size) {
    int end_of_file_flag = 0;
    do {
      if (!std::fgets(next_line_pointer, max_line_size, file_pointer)) {
        end_of_file_flag = 1;
        break;
      }

      while (next_line_pointer[0] == ' ' || next_line_pointer[0] == '\t' || next_line_pointer[0] == '\n'
             || next_line_pointer[0] == '\r' || next_line_pointer[0] == '\f') {
        next_line_pointer++;
      }

    } while ((std::strncmp("#", next_line_pointer, 1) == 0) || (std::strlen(next_line_pointer) == 0));

    // remove comments starting with `#' in a line
    char *pch = std::strchr(next_line_pointer, '#');
    if (pch) {
      *pch = '\0';
    }

    return end_of_file_flag;
  }

  inline int Utility::GetNumberOfWordsInLine(char const *line_pointer) {
    int numWords = 0;

    char const *buf = line_pointer;
    char c = *buf;
    while (c) {
      if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
        c = *++buf;
        continue;
      };

      ++numWords;
      c = *++buf;

      while (c) {
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
          break;
        }
        c = *++buf;
      }
    }

    return numWords;
  }

  inline std::vector<std::string> Utility::GetWordsInLine(char const *line_pointer) {
    char const *buf = line_pointer;
    std::string line_pointer_string(buf);
    std::vector<std::string> words;
    std::size_t beg = 0;
    std::size_t len = 0;
    std::size_t add = 0;
    char c = *buf;

    while (c) {
      // leading whitespace
      if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
        c = *++buf;
        ++beg;
        continue;
      };

      len = 0;

    // handle escaped/quoted text.
    quoted:

      // handle single quote
      if (c == '\'') {
        ++beg;
        add = 1;
        c = *++buf;
        while (((c != '\'') && (c != '\0')) || ((c == '\\') && (buf[1] == '\''))) {
          if ((c == '\\') && (buf[1] == '\'')) {
            ++buf;
            ++len;
          }
          c = *++buf;
          ++len;
        }
        if (c != '\'') {
          ++len;
        }
        c = *++buf;

        // handle double quote
      } else if (c == '"') {
        ++beg;
        add = 1;
        c = *++buf;
        while (((c != '"') && (c != '\0')) || ((c == '\\') && (buf[1] == '"'))) {
          if ((c == '\\') && (buf[1] == '"')) {
            ++buf;
            ++len;
          }
          c = *++buf;
          ++len;
        }
        if (c != '"') {
          ++len;
        }
        c = *++buf;
      }

      // unquoted
      while (1) {
        if ((c == '\'') || (c == '"')) {
          goto quoted;
        }
        // skip escaped quote
        if ((c == '\\') && ((buf[1] == '\'') || (buf[1] == '"'))) {
          ++buf;
          ++len;
          c = *++buf;
          ++len;
        }
        if ((c == ' ') || (c == '\t') || (c == '\r') || (c == '\n') || (c == '\f') || (c == '\0')) {
          words.push_back(line_pointer_string.substr(beg, len));
          beg += len + add;
          break;
        }
        c = *++buf;
        ++len;
      }
    }
    return words;
  }

}  // namespace openKIM

#endif  // UTILITY_HPP