//
// helper.hpp
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
// Copyright (c) 2019--2021, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Yaser Afshar
//

#ifndef HELPER_HPP
#define HELPER_HPP

#pragma GCC system_header

/*!
 * \file helper.hpp
 *
 * \author Yaser Afshar (yafshar@umn.edu)
 *
 * \brief This file contains helper classes for memory management
 *        which are STL like extension for multi-dimensional arrays
 *
 */

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace openKIM {

#ifdef HELPER_LOG_ERROR
#  undef HELPER_LOG_ERROR
#endif

#ifdef HELPER_LOG_WARNING
#  undef HELPER_LOG_WARNING
#endif

/*!
 * \brief Helper macro for printing error message
 * formats messages, filename, line number and function
 * name into an std::ostringstream object
 *
 */
#define HELPER_LOG_ERROR(msg)                                                                            \
  {                                                                                                      \
    std::ostringstream ss;                                                                               \
    ss << "\nError :" << __FILE__ << ":" << __LINE__ << ":@(" << __FUNCTION__ << ")\n" << msg << "\n\n"; \
    std::cerr << ss.str();                                                                               \
  }

#define HELPER_LOG_WARNING(msg)                                                                            \
  {                                                                                                        \
    std::ostringstream ss;                                                                                 \
    ss << "\nWarning :" << __FILE__ << ":" << __LINE__ << ":@(" << __FUNCTION__ << ")\n" << msg << "\n\n"; \
    std::cerr << ss.str();                                                                                 \
  }

  /*! Type alias for vector of constant dimension 3 */
  using VectorOfSizeDIM = double[3];

  /*! Type alias for vector of constant dimension 6 */
  using VectorOfSizeSix = double[6];

  /*! \class _ArrayBasic The basic STL like container similar to `std::vector` to
   * handle multi-dimensional arrays.
   *
   * \brief An STL like container similar to <a
   * href="https://en.cppreference.com/w/cpp/container/vector">std::vector</a>
   * that encapsulates dynamic size arrays in a sequence container
   *
   * \tparam DataType The type of the element. Default (double)
   */
  template <class DataType>
  class _ArrayBasic {
  public:
    /*!
     * \brief Construct a new _ArrayBasic object
     *
     */
    _ArrayBasic();

    /*!
     * \brief Construct a new _ArrayBasic object
     *
     * \param count The size of the container
     */
    explicit _ArrayBasic(std::size_t const count);

    /*!
     * \brief Construct a new _ArrayBasic object
     *
     * \param count The size of the container
     * \param value The value to initialize element of the container with
     */
    _ArrayBasic(std::size_t const count, DataType const value);

    /*!
     * \brief Construct a new _ArrayBasic object
     *
     * \param count The size of the container
     * \param array The array of data to initialize element of the container
     * with.
     */
    _ArrayBasic(std::size_t const count, DataType const *array);

    /*!
     * \brief Construct a new _ArrayBasic object
     * Copy constructor. Constructs the container with the copy of the contents of
     * other.
     *
     * \param other Another container to be used as source to initialize the
     * element of the container with
     */
    _ArrayBasic(_ArrayBasic<DataType> const &other);

    /*!
     * \brief Construct a new _ArrayBasic object
     * Move constructor. Constructs the container with the contents of other using
     * move semantics.
     *
     * \param other Another container to be used as source to initialize the
     * element of the container with
     */
    _ArrayBasic(_ArrayBasic<DataType> &&other);

    /*!
     * \brief Destroy the _ArrayBasic object
     *
     */
    ~_ArrayBasic();

    /*!
     * \brief Copy assignment operator. Replaces the contents with a copy of the
     * contents of other
     *
     * \param other Another container to use as data source
     *
     * \return _ArrayBasic<DataType>&
     */
    _ArrayBasic<DataType> &operator=(_ArrayBasic<DataType> const &other);

    /*!
     * \brief Move assignment operator. Replaces the contents with those of other
     * using move semantics
     *
     * \param other Another container to use as data source
     *
     * \return _ArrayBasic<DataType>&
     */
    _ArrayBasic<DataType> &operator=(_ArrayBasic<DataType> &&other);

    /*!
     * \brief Returns pointer to the underlying array serving as element storage.
     *
     * \return DataType*
     */
    inline DataType const *data() const noexcept;
    inline DataType *data() noexcept;

    /*!
     * \brief Returns the number of element in the container.
     *
     * \return std::size_t
     */
    inline std::size_t size() const;

    /*!
     * \brief Erases all element from the container.
     *
     */
    inline void clear() noexcept;

    /*!
     * \brief Requests the removal of unused capacity.
     *
     */
    inline void shrink_to_fit();

    /*!
     * \brief Returns the number of element that the container has currently
     * allocated space for.
     *
     * \return std::size_t
     */
    inline std::size_t capacity() const noexcept;

    /*!
     * \brief Appends the given element value to the end of the container.
     *
     * \param value
     */
    inline void push_back(DataType const &value);
    inline void push_back(DataType &&value);

  protected:
    /*!
     * \brief Check the index range based on container size
     *
     * \param n Index
     */
    inline void _range_check(int n) const;

    /*!
     * \brief Check the index range based on input size
     *
     * \param n Index
     * \param tsize Input size
     */
    inline void _range_check(int n, std::size_t tsize) const;

  protected:
    /*! Dynamic contiguous array */
    std::vector<DataType> m_;
  };

  /*! \class Array1DView A 1-dimensional STL like container.
   *
   * \brief A 1-dimensional STL like container that encapsulates dynamic size
   * arrays in a sequence container
   *
   * \tparam DataType The type of the element. Default (double)
   */
  template <class DataType>
  class Array1DView {
  public:
    /*!
     * \brief Construct a new Array1DView object
     *
     * \param count The size of the container
     * \param array The array of data to initialize element of the container
     * with.
     */
    Array1DView(std::size_t const count, DataType *array);
    Array1DView(std::size_t const count, DataType const *array);

    /*!
     * \brief Construct a new Array1DView object
     * Copy constructor. Constructs the container with the copy of the contents of
     * other.
     *
     * \param other Another container to be used as source to initialize the
     * element of the container with
     */
    Array1DView(Array1DView<DataType> const &other);

    /*!
     * \brief Destroy the Array1DView object
     *
     */
    ~Array1DView();

    /*!
     * \brief Returns pointer to the underlying array serving as element storage.
     *
     * \return DataType*
     */
    inline DataType const *data() const noexcept;
    inline DataType *data() noexcept;

    /*!
     * \brief Returns the element at specified location \b i.
     * No bounds checking is performed.
     *
     * \param i Position of the element to return
     *
     * \return const DataType The requested element.
     */
    inline const DataType operator()(int i) const;
    inline DataType &operator()(int i);

    /*!
     * \brief Returns the element at specified location \b i , with bounds
     * checking.
     *
     * \param i Position of the element to return
     *
     * \return const DataType The requested element.
     */
    inline DataType const at(int i) const;
    inline DataType &at(int i);

    /*!
     * \brief Returns the element at specified location \b i.
     * No bounds checking is performed.
     *
     * \param i Position of the element to return
     *
     * \return const DataType The requested element.
     */
    inline const DataType operator[](int i) const;
    inline DataType &operator[](int i);

  private:
    Array1DView() = delete;

    Array1DView<DataType> &operator=(Array1DView<DataType> const &other) = delete;

    Array1DView<DataType> &operator=(Array1DView<DataType> &&other) = delete;

  protected:
    /*!
     * \brief Check the index range based on input size
     *
     * \param n Index
     * \param tsize Input size
     */
    inline void _range_check(int n, std::size_t tsize) const;

  protected:
    /*! The extent of the container in the 1st mode */
    std::size_t extent_zero_;

    /*! Data pointer */
    DataType *const m_;
  };

  template <class DataType>
  class Array2DView {
  public:
    /*!
     * \brief Construct a new Array2DView object
     *
     * \param extent_zero The extent of the container in the 1st mode
     * \param extent_one The extent of the container in the 2nd mode
     * \param array The array of data to set the pointer of the container to it.
     */
    Array2DView(std::size_t const extent_zero, std::size_t const extent_one, DataType *array);

    Array2DView(std::size_t const extent_zero, std::size_t const extent_one, DataType const *array);

    /*!
     * \brief Construct a new Array2DView object
     * Copy constructor. Constructs the container with the copy of the contents of
     * other.
     *
     * \param other Another container to be used as source to initialize the
     * element of the container with
     */
    Array2DView(Array2DView<DataType> const &other);

    /*!
     * \brief Destroy the Array2D object
     *
     */
    ~Array2DView();

    /*!
     * \brief Returns pointer to the underlying array serving as element storage.
     *
     * \return DataType*
     */
    inline DataType const *data() const noexcept;
    inline DataType *data() noexcept;

    inline Array1DView<DataType> const data_1D(int i) const;
    inline Array1DView<DataType> data_1D(int i);

    /*!
     * \brief Returns the element at specified location \b (i, j).
     * No bounds checking is performed.
     *
     * \param i Position of the element in the 1st mode
     * \param j Position of the element in the 2nd mode
     *
     * \return const DataType The requested element.
     */
    inline const DataType operator()(int i, int j) const;
    inline DataType &operator()(int i, int j);

    /*!
     * \brief Returns the element at specified location \b (i, j) , with bounds
     * checking.
     *
     * \param i Position of the element in the 1st mode
     * \param j Position of the element in the 2nd mode
     *
     * \return const DataType The requested element.
     */
    inline DataType const at(int i, int j) const;
    inline DataType &at(int i, int j);

    /*! \class j_operator A helper class to provide multidimensional array access
     * semantics.
     *
     * \brief To provide 2-dimensional array access semantics, operator[] has to
     * return a reference to a 1D vector, which has to have its own operator[]
     * which returns a reference to the element.
     */
    class j_operator {
    public:
      /*!
       * \brief Construct a new j_operator object
       *
       * \param array Refernce to Array2D class
       * \param i Position of the element in the 1st mode
       */
      j_operator(Array2DView<DataType> &array, int i);

      /*!
       * \brief Provide array-like access and returns the element at specified
       * location \b [i][j]. No bounds checking is performed.
       *
       * \param j Position of the element in the 2nd mode
       *
       * \return const DataType The requested element.
       */
      inline const DataType operator[](int j) const;
      inline DataType &operator[](int j);

    private:
      /*! Refernce to Array2D class */
      Array2DView<DataType> &j_array_;

      std::size_t i_;
    };

    /*!
     * \brief Provide array-like access and returns the element at specified
     * location \b [i][j]. No bounds checking is performed.
     *
     * \param i Position of the element in the 1st mode
     * \param j Position of the element in the 2nd mode
     *
     * \return const DataType The requested element.
     *
     * \note
     * To provide multidimensional array access semantics, we are using multiple
     * overloads for \code operator[] \endcode . For speed kOne should avoid this
     * complexity, uses \code operator() \endcode as \code (i, j) \endcode
     * directly.
     */
    inline const j_operator operator[](int i) const;
    inline j_operator operator[](int i);

  private:
    Array2DView() = delete;

    Array2DView<DataType> &operator=(Array2DView<DataType> const &other) = delete;

    Array2DView<DataType> &operator=(Array2DView<DataType> &&other) = delete;

  protected:
    /*!
     * \brief Check the index range based on input size
     *
     * \param n Index
     * \param tsize Input size
     */
    inline void _range_check(int n, std::size_t tsize) const;

  protected:
    /*! The extent of the container in the 1st mode */
    std::size_t extent_zero_;

    /*! The extent of the container in the 2nd mode */
    std::size_t extent_one_;

    /*! Data pointer */
    DataType *const m_;
  };

  /*! \class Array2D A 2-dimensional STL like container.
   *
   * \brief A 2-dimensional STL like container that encapsulates dynamic size
   * arrays with a 2-dimensional shape in a sequence container
   *
   * \tparam DataType The type of the element. Default (double)
   */
  template <class DataType>
  class Array2D : public _ArrayBasic<DataType> {
  public:
    /*!
     * \brief Construct a new Array2D object
     *
     */
    Array2D();

    /*!
     * \brief Construct a new Array2D object
     *
     * \param extent_zero The extent of the container in the 1st mode
     * \param extent_one The extent of the container in the 2nd mode
     */
    Array2D(std::size_t const extent_zero, std::size_t const extent_one);

    /*!
     * \brief Construct a new Array2D object
     *
     * \param extent_zero The extent of the container in the 1st mode
     * \param extent_one The extent of the container in the 2nd mode
     * \param value The value to initialize element of the container with
     */
    Array2D(std::size_t const extent_zero, std::size_t const extent_one, DataType const value);

    /*!
     * \brief Construct a new Array2D object
     *
     * \param extent_zero The extent of the container in the 1st mode
     * \param extent_one The extent of the container in the 2nd mode
     * \param array The array of data to initialize element of the container in a
     * row-major format.
     */
    Array2D(std::size_t const extent_zero, std::size_t const extent_one, DataType const *array);

    /*!
     * \brief Construct a new Array2D object
     * Copy constructor. Constructs the container with the copy of the contents of
     * other.
     *
     * \param other Another container to be used as source to initialize the
     * element of the container with
     */
    Array2D(Array2D<DataType> const &other);

    /*!
     * \brief Construct a new Array2D object
     * Move constructor. Constructs the container with the contents of other using
     * move semantics.
     *
     * \param other Another container to be used as source to initialize the
     * element of the container with
     */
    Array2D(Array2D<DataType> &&other);

    /*!
     * \brief Destroy the Array2D object
     *
     */
    ~Array2D();

    /*!
     * \brief Copy assignment operator. Replaces the contents with a copy of the
     * contents of other
     *
     * \param other Another container to use as data source
     *
     * \return Array2D<DataType>&
     */
    Array2D<DataType> &operator=(Array2D<DataType> const &other);

    /*!
     * \brief Move assignment operator. Replaces the contents with those of other
     * using move semantics
     *
     * \param other Another container to use as data source
     *
     * \return Array2D<DataType>&
     */
    Array2D<DataType> &operator=(Array2D<DataType> &&other);

    /*!
     * \brief Returns Array1DView to the underlying starting element at row \c i.
     * \sa Array1DView
     *
     * \param i Row index
     * \return Array1DView<DataType>
     */
    inline Array1DView<DataType> const data_1D(int i) const;
    inline Array1DView<DataType> data_1D(int i);

    /*!
     * \brief Resizes the container to contain \c extent_zero times \c extent_one
     * element.
     *
     * \param extent_zero
     * \param extent_one
     */
    inline void resize(int const extent_zero, int const extent_one);

    /*!
     * \brief Resizes the container to contain \c extent_zero times \c extent_one
     * element.
     *
     * \param extent_zero
     * \param extent_one
     * \param new_value The new value to initialize the new element with
     */
    inline void resize(int const extent_zero, int const extent_one, DataType const new_value);

    /*!
     * \brief Resizes the container to contain \c extent_zero times \c extent_one
     * element.
     *
     * \param extent_zero
     * \param extent_one
     * \param new_array The new array of data to initialize element of the
     * container with.
     */
    inline void resize(int const extent_zero, int const extent_one, DataType const *new_array);

    /*!
     * \brief Returns the element at specified location \b (i, j).
     * No bounds checking is performed.
     *
     * \param i Position of the element in the 1st mode
     * \param j Position of the element in the 2nd mode
     *
     * \return const DataType The requested element.
     */
    inline const DataType operator()(int i, int j) const;
    inline DataType &operator()(int i, int j);

    /*!
     * \brief Returns the element at specified location \b (i, j) , with bounds
     * checking.
     *
     * \param i Position of the element in the 1st mode
     * \param j Position of the element in the 2nd mode
     *
     * \return const DataType The requested element.
     */
    inline DataType const at(int i, int j) const;
    inline DataType &at(int i, int j);

    /*! \class j_operator A helper class to provide multidimensional array access
     * semantics.
     *
     * \brief To provide 2-dimensional array access semantics, operator[] has to
     * return a reference to a 1D vector, which has to have its own operator[]
     * which returns a reference to the element.
     */
    class j_operator {
    public:
      /*!
       * \brief Construct a new j_operator object
       *
       * \param array Refernce to Array2D class
       * \param i Position of the element in the 1st mode
       */
      j_operator(Array2D<DataType> &array, int i);

      /*!
       * \brief Provide array-like access and returns the element at specified
       * location \b [i][j]. No bounds checking is performed.
       *
       * \param j Position of the element in the 2nd mode
       *
       * \return const DataType The requested element.
       */
      inline const DataType operator[](int j) const;
      inline DataType &operator[](int j);

    private:
      /*! Refernce to Array2D class */
      Array2D<DataType> &j_array_;

      std::size_t i_;
    };

    /*!
     * \brief Provide array-like access and returns the element at specified
     * location \b [i][j]. No bounds checking is performed.
     *
     * \param i Position of the element in the 1st mode
     * \param j Position of the element in the 2nd mode
     *
     * \return const DataType The requested element.
     *
     * \note
     * To provide multidimensional array access semantics, we are using multiple
     * overloads for \code operator[] \endcode . For speed kOne should avoid this
     * complexity, uses \code operator() \endcode as \code (i, j) \endcode
     * directly.
     */
    inline const j_operator operator[](int i) const;
    inline j_operator operator[](int i);

  protected:
    /*! The extent of the container in the 1st mode */
    std::size_t extent_zero_;

    /*! The extent of the container in the 2nd mode */
    std::size_t extent_one_;
  };

  template <class DataType>
  class Array3D : public _ArrayBasic<DataType> {
  public:
    Array3D();

    Array3D(std::size_t const extent_zero, std::size_t const extent_one, std::size_t const extent_two);

    Array3D(std::size_t const extent_zero,
            std::size_t const extent_one,
            std::size_t const extent_two,
            DataType const value);

    Array3D(std::size_t const extent_zero,
            std::size_t const extent_one,
            std::size_t const extent_two,
            DataType const *array);

    Array3D(Array3D<DataType> const &other);

    Array3D(Array3D<DataType> &&other);

    ~Array3D();

    Array3D<DataType> &operator=(Array3D<DataType> const &other);

    Array3D<DataType> &operator=(Array3D<DataType> &&other);

    inline Array2DView<DataType> const data_2D(int i) const;
    inline Array2DView<DataType> data_2D(int i);

    inline Array1DView<DataType> const data_1D(int i, int j) const;
    inline Array1DView<DataType> data_1D(int i, int j);

    inline void resize(int const extent_zero, int const extent_one, int const extent_two);

    inline void resize(int const extent_zero,
                       int const extent_one,
                       int const extent_two,
                       DataType const new_value);

    inline void resize(int const extent_zero,
                       int const extent_one,
                       int const extent_two,
                       DataType const *new_array);

    inline const DataType operator()(int i, int j, int k) const;
    inline DataType &operator()(int i, int j, int k);

    inline DataType const at(int i, int j, int k) const;
    inline DataType &at(int i, int j, int k);

    class j_operator {
    public:
      j_operator(Array3D<DataType> &array, int i);

      class k_operator {
      public:
        k_operator(Array3D<DataType> &array, int i, int j);

        inline const DataType operator[](int k) const;
        inline DataType &operator[](int k);

      private:
        Array3D<DataType> &k_array_;

        std::size_t i_;
        std::size_t j_;
      };

      inline const k_operator operator[](int j) const;
      inline k_operator operator[](int j);

    private:
      Array3D<DataType> &j_array_;

      std::size_t i_;
    };

    inline const j_operator operator[](int i) const;
    inline j_operator operator[](int i);

  protected:
    std::size_t extent_zero_;
    std::size_t extent_one_;
    std::size_t extent_two_;
  };

  template <class DataType>
  _ArrayBasic<DataType>::_ArrayBasic() {}

  template <class DataType>
  _ArrayBasic<DataType>::_ArrayBasic(std::size_t const count) : m_(count) {}

  template <class DataType>
  _ArrayBasic<DataType>::_ArrayBasic(std::size_t const count, DataType const value) : m_(count, value) {}

  template <class DataType>
  _ArrayBasic<DataType>::_ArrayBasic(std::size_t const count, DataType const *array)
      : m_(array, array + count) {}

  template <class DataType>
  _ArrayBasic<DataType>::_ArrayBasic(_ArrayBasic<DataType> const &other) : m_(other.m_) {}

  template <class DataType>
  _ArrayBasic<DataType>::_ArrayBasic(_ArrayBasic<DataType> &&other) : m_(std::move(other.m_)) {}

  template <class DataType>
  _ArrayBasic<DataType>::~_ArrayBasic() {}

  template <class DataType>
  _ArrayBasic<DataType> &_ArrayBasic<DataType>::operator=(_ArrayBasic<DataType> const &other) {
    m_.resize(other.size());
    std::copy(other.m_.begin(), other.m_.end(), m_.begin());
    return *this;
  }

  template <class DataType>
  _ArrayBasic<DataType> &_ArrayBasic<DataType>::operator=(_ArrayBasic<DataType> &&other) {
    m_ = std::move(other.m_);
    return *this;
  }

  template <class DataType>
  inline DataType const *_ArrayBasic<DataType>::data() const noexcept {
    return m_.data();
  }

  template <class DataType>
  inline DataType *_ArrayBasic<DataType>::data() noexcept {
    return m_.data();
  }

  template <class DataType>
  inline std::size_t _ArrayBasic<DataType>::size() const {
    return m_.size();
  }

  template <class DataType>
  inline void _ArrayBasic<DataType>::clear() noexcept {
    m_.clear();
  }

  template <class DataType>
  inline void _ArrayBasic<DataType>::shrink_to_fit() {
    m_.shrink_to_fit();
  }

  template <class DataType>
  inline std::size_t _ArrayBasic<DataType>::capacity() const noexcept {
    return m_.capacity();
  }

  template <class DataType>
  inline void _ArrayBasic<DataType>::push_back(DataType const &value) {
    m_.push_back(value);
  }

  template <class DataType>
  inline void _ArrayBasic<DataType>::push_back(DataType &&value) {
    m_.push_back(value);
  }

  template <class DataType>
  inline void _ArrayBasic<DataType>::_range_check(int n) const {
    if (n >= size()) {
      HELPER_LOG_ERROR("The input index is out of range! " + std::to_string(n)
                       + " >= " + std::to_string(size()));
      std::abort();
    }
  }

  template <class DataType>
  inline void _ArrayBasic<DataType>::_range_check(int n, std::size_t tsize) const {
    if (n >= tsize) {
      HELPER_LOG_ERROR("The input index is out of range! " + std::to_string(n)
                       + " >= " + std::to_string(tsize));
      std::abort();
    }
  }

  template <class DataType>
  Array1DView<DataType>::Array1DView(std::size_t const count, DataType *array)
      : extent_zero_(count), m_(array) {}

  template <class DataType>
  Array1DView<DataType>::Array1DView(std::size_t const count, DataType const *array)
      : extent_zero_(count), m_(const_cast<DataType *>(array)) {}

  template <class DataType>
  Array1DView<DataType>::Array1DView(Array1DView<DataType> const &other)
      : extent_zero_(other.extent_zero_), m_(other.m_) {}

  template <class DataType>
  Array1DView<DataType>::~Array1DView() {}

  template <class DataType>
  inline DataType const *Array1DView<DataType>::data() const noexcept {
    return m_;
  }

  template <class DataType>
  inline DataType *Array1DView<DataType>::data() noexcept {
    return m_;
  }

  template <class DataType>
  inline const DataType Array1DView<DataType>::operator()(int i) const {
    return m_[i];
  }

  template <class DataType>
  inline DataType &Array1DView<DataType>::operator()(int i) {
    return m_[i];
  }

  template <class DataType>
  inline DataType &Array1DView<DataType>::at(int i) {
    _range_check(i, extent_zero_);
    return m_[i];
  }

  template <class DataType>
  inline DataType const Array1DView<DataType>::at(int i) const {
    _range_check(i, extent_zero_);
    return m_[i];
  }

  template <class DataType>
  inline const DataType Array1DView<DataType>::operator[](int i) const {
    return m_[i];
  }

  template <class DataType>
  inline DataType &Array1DView<DataType>::operator[](int i) {
    return m_[i];
  }

  template <class DataType>
  inline void Array1DView<DataType>::_range_check(int n, std::size_t tsize) const {
    if (n >= tsize) {
      HELPER_LOG_ERROR("The input index is out of range! " + std::to_string(n)
                       + " >= " + std::to_string(tsize));
      std::abort();
    }
  }

  template <class DataType>
  Array2DView<DataType>::Array2DView(std::size_t const extent_zero,
                                     std::size_t const extent_one,
                                     DataType *array)
      : extent_zero_(extent_zero), extent_one_(extent_one), m_(array) {}

  template <class DataType>
  Array2DView<DataType>::Array2DView(std::size_t const extent_zero,
                                     std::size_t const extent_one,
                                     DataType const *array)
      : extent_zero_(extent_zero), extent_one_(extent_one), m_(const_cast<DataType *>(array)) {}

  template <class DataType>
  Array2DView<DataType>::Array2DView(Array2DView<DataType> const &other)
      : extent_zero_(other.extent_zero_), extent_one_(other.extent_one_), m_(other.m_) {}

  template <class DataType>
  Array2DView<DataType>::~Array2DView() {}

  template <class DataType>
  inline DataType const *Array2DView<DataType>::data() const noexcept {
    return m_;
  }

  template <class DataType>
  inline DataType *Array2DView<DataType>::data() noexcept {
    return m_;
  }

  template <class DataType>
  inline Array1DView<DataType> const Array2DView<DataType>::data_1D(int i) const {
    return Array1DView<DataType>(extent_one_, m_ + i * extent_one_);
  }

  template <class DataType>
  inline Array1DView<DataType> Array2DView<DataType>::data_1D(int i) {
    return Array1DView<DataType>(extent_one_, m_ + i * extent_one_);
  }

  template <class DataType>
  inline const DataType Array2DView<DataType>::operator()(int i, int j) const {
    std::size_t const n = i * extent_one_ + j;
    return m_[n];
  }

  template <class DataType>
  inline DataType &Array2DView<DataType>::operator()(int i, int j) {
    std::size_t const n = i * extent_one_ + j;
    return m_[n];
  }

  template <class DataType>
  inline DataType &Array2DView<DataType>::at(int i, int j) {
    _range_check(i, extent_zero_);
    _range_check(j, extent_one_);
    std::size_t const n = i * extent_one_ + j;
    return m_[n];
  }

  template <class DataType>
  inline DataType const Array2DView<DataType>::at(int i, int j) const {
    _range_check(i, extent_zero_);
    _range_check(j, extent_one_);
    std::size_t const n = i * extent_one_ + j;
    return m_[n];
  }

  template <class DataType>
  Array2DView<DataType>::j_operator::j_operator(Array2DView<DataType> &array, int i)
      : j_array_(array), i_(i) {}

  template <class DataType>
  inline const DataType Array2DView<DataType>::j_operator::operator[](int j) const {
    std::size_t const n = i_ * j_array_.extent_one_ + j;
    return j_array_.m_[n];
  }

  template <class DataType>
  inline DataType &Array2DView<DataType>::j_operator::operator[](int j) {
    std::size_t const n = i_ * j_array_.extent_one_ + j;
    return j_array_.m_[n];
  }

  template <class DataType>
  inline const typename Array2DView<DataType>::j_operator Array2DView<DataType>::operator[](int i) const {
    return j_operator(*this, i);
  }

  template <class DataType>
  inline typename Array2DView<DataType>::j_operator Array2DView<DataType>::operator[](int i) {
    return j_operator(*this, i);
  }

  template <class DataType>
  inline void Array2DView<DataType>::_range_check(int n, std::size_t tsize) const {
    if (n >= tsize) {
      HELPER_LOG_ERROR("The input index is out of range! " + std::to_string(n)
                       + " >= " + std::to_string(tsize));
      std::abort();
    }
  }

  template <class DataType>
  Array2D<DataType>::Array2D() : _ArrayBasic<DataType>(), extent_zero_(0), extent_one_(0) {}

  template <class DataType>
  Array2D<DataType>::Array2D(std::size_t const extent_zero, std::size_t const extent_one)
      : _ArrayBasic<DataType>(extent_zero * extent_one), extent_zero_(extent_zero), extent_one_(extent_one) {}

  template <class DataType>
  Array2D<DataType>::Array2D(std::size_t const extent_zero,
                             std::size_t const extent_one,
                             DataType const value)
      : _ArrayBasic<DataType>(extent_zero * extent_one, value),
        extent_zero_(extent_zero),
        extent_one_(extent_one) {}

  template <class DataType>
  Array2D<DataType>::Array2D(std::size_t const extent_zero,
                             std::size_t const extent_one,
                             DataType const *array)
      : _ArrayBasic<DataType>(extent_zero * extent_one, array),
        extent_zero_(extent_zero),
        extent_one_(extent_one) {}

  template <class DataType>
  Array2D<DataType>::Array2D(Array2D<DataType> const &other)
      : _ArrayBasic<DataType>(other), extent_zero_(other.extent_zero_), extent_one_(other.extent_one_) {}

  template <class DataType>
  Array2D<DataType>::Array2D(Array2D<DataType> &&other)
      : _ArrayBasic<DataType>(std::move(other)),
        extent_zero_(other.extent_zero_),
        extent_one_(other.extent_one_) {}

  template <class DataType>
  Array2D<DataType>::~Array2D() {}

  template <class DataType>
  Array2D<DataType> &Array2D<DataType>::operator=(Array2D<DataType> const &other) {
    _ArrayBasic<DataType>::operator=(other);
    extent_zero_ = other.extent_zero_;
    extent_one_ = other.extent_one_;
    return *this;
  }

  template <class DataType>
  Array2D<DataType> &Array2D<DataType>::operator=(Array2D<DataType> &&other) {
    _ArrayBasic<DataType>::operator=(std::move(other));
    extent_zero_ = other.extent_zero_;
    extent_one_ = other.extent_one_;
    return *this;
  }

  template <class DataType>
  inline Array1DView<DataType> const Array2D<DataType>::data_1D(int i) const {
    return Array1DView<DataType>(extent_one_, this->m_.data() + i * extent_one_);
  }

  template <class DataType>
  inline Array1DView<DataType> Array2D<DataType>::data_1D(int i) {
    return Array1DView<DataType>(extent_one_, this->m_.data() + i * extent_one_);
  }

  template <class DataType>
  inline void Array2D<DataType>::resize(int const extent_zero, int const extent_one) {
    extent_zero_ = extent_zero;
    extent_one_ = extent_one;
    std::size_t const n = extent_zero_ * extent_one_;
    this->m_.resize(n);
  }

  template <class DataType>
  inline void Array2D<DataType>::resize(int const extent_zero,
                                        int const extent_one,
                                        DataType const new_value) {
    extent_zero_ = extent_zero;
    extent_one_ = extent_one;
    std::size_t const n = extent_zero_ * extent_one_;
    this->m_.resize(n, new_value);
  }

  template <class DataType>
  inline void Array2D<DataType>::resize(int const extent_zero,
                                        int const extent_one,
                                        DataType const *new_array) {
    extent_zero_ = extent_zero;
    extent_one_ = extent_one;
    std::size_t const n = extent_zero_ * extent_one_;
    this->m_.resize(n);
    std::copy(new_array, new_array + n, this->m_.data());
  }

  template <class DataType>
  inline const DataType Array2D<DataType>::operator()(int i, int j) const {
    std::size_t const n = i * extent_one_ + j;
    return this->m_[n];
  }

  template <class DataType>
  inline DataType &Array2D<DataType>::operator()(int i, int j) {
    std::size_t const n = i * extent_one_ + j;
    return this->m_[n];
  }

  template <class DataType>
  inline DataType &Array2D<DataType>::at(int i, int j) {
    this->_range_check(i, extent_zero_);
    this->_range_check(j, extent_one_);
    std::size_t const n = i * extent_one_ + j;
    return this->m_[n];
  }

  template <class DataType>
  inline DataType const Array2D<DataType>::at(int i, int j) const {
    this->_range_check(i, extent_zero_);
    this->_range_check(j, extent_one_);
    std::size_t const n = i * extent_one_ + j;
    return this->m_[n];
  }

  template <class DataType>
  Array2D<DataType>::j_operator::j_operator(Array2D<DataType> &array, int i) : j_array_(array), i_(i) {}

  template <class DataType>
  inline const DataType Array2D<DataType>::j_operator::operator[](int j) const {
    std::size_t const n = i_ * j_array_.extent_one_ + j;
    return j_array_.m_[n];
  }

  template <class DataType>
  inline DataType &Array2D<DataType>::j_operator::operator[](int j) {
    std::size_t const n = i_ * j_array_.extent_one_ + j;
    return j_array_.m_[n];
  }

  template <class DataType>
  inline const typename Array2D<DataType>::j_operator Array2D<DataType>::operator[](int i) const {
    return j_operator(*this, i);
  }

  template <class DataType>
  inline typename Array2D<DataType>::j_operator Array2D<DataType>::operator[](int i) {
    return j_operator(*this, i);
  }

  template <class DataType>
  Array3D<DataType>::Array3D() : _ArrayBasic<DataType>(), extent_zero_(0), extent_one_(0), extent_two_(0) {}

  template <class DataType>
  Array3D<DataType>::Array3D(std::size_t const extent_zero,
                             std::size_t const extent_one,
                             std::size_t const extent_two)
      : _ArrayBasic<DataType>(extent_zero * extent_one * extent_two),
        extent_zero_(extent_zero),
        extent_one_(extent_one),
        extent_two_(extent_two) {}

  template <class DataType>
  Array3D<DataType>::Array3D(std::size_t const extent_zero,
                             std::size_t const extent_one,
                             std::size_t const extent_two,
                             DataType const value)
      : _ArrayBasic<DataType>(extent_zero * extent_one * extent_two, value),
        extent_zero_(extent_zero),
        extent_one_(extent_one),
        extent_two_(extent_two) {}

  template <class DataType>
  Array3D<DataType>::Array3D(std::size_t const extent_zero,
                             std::size_t const extent_one,
                             std::size_t const extent_two,
                             DataType const *array)
      : _ArrayBasic<DataType>(extent_zero * extent_one * extent_two, array),
        extent_zero_(extent_zero),
        extent_one_(extent_one),
        extent_two_(extent_two) {}

  template <class DataType>
  Array3D<DataType>::Array3D(Array3D<DataType> const &other)
      : _ArrayBasic<DataType>(other),
        extent_zero_(other.extent_zero_),
        extent_one_(other.extent_one_),
        extent_two_(other.extent_two_) {}

  template <class DataType>
  Array3D<DataType>::Array3D(Array3D<DataType> &&other)
      : _ArrayBasic<DataType>(std::move(other)),
        extent_zero_(other.extent_zero_),
        extent_one_(other.extent_one_),
        extent_two_(other.extent_two_) {}

  template <class DataType>
  Array3D<DataType>::~Array3D() {}

  template <class DataType>
  Array3D<DataType> &Array3D<DataType>::operator=(Array3D<DataType> const &other) {
    _ArrayBasic<DataType>::operator=(other);
    extent_zero_ = other.extent_zero_;
    extent_one_ = other.extent_one_;
    extent_two_ = other.extent_two_;
    return *this;
  }

  template <class DataType>
  Array3D<DataType> &Array3D<DataType>::operator=(Array3D<DataType> &&other) {
    _ArrayBasic<DataType>::operator=(std::move(other));
    extent_zero_ = other.extent_zero_;
    extent_one_ = other.extent_one_;
    extent_two_ = other.extent_two_;
    return *this;
  }

  template <class DataType>
  inline Array2DView<DataType> const Array3D<DataType>::data_2D(int i) const {
    return Array2DView<DataType>(extent_one_, extent_two_, this->m_.data() + i * extent_one_ * extent_two_);
  }

  template <class DataType>
  inline Array2DView<DataType> Array3D<DataType>::data_2D(int i) {
    return Array2DView<DataType>(extent_one_, extent_two_, this->m_.data() + i * extent_one_ * extent_two_);
  }

  template <class DataType>
  inline Array1DView<DataType> const Array3D<DataType>::data_1D(int i, int j) const {
    return Array1DView<DataType>(extent_two_, this->m_.data() + (i * extent_one_ + j) * extent_two_);
  }

  template <class DataType>
  inline Array1DView<DataType> Array3D<DataType>::data_1D(int i, int j) {
    return Array1DView<DataType>(extent_two_, this->m_.data() + (i * extent_one_ + j) * extent_two_);
  }

  template <class DataType>
  inline void Array3D<DataType>::resize(int const extent_zero, int const extent_one, int const extent_two) {
    extent_zero_ = extent_zero;
    extent_one_ = extent_one;
    extent_two_ = extent_two;
    std::size_t const n = extent_zero_ * extent_one_ * extent_two_;
    this->m_.resize(n);
  }

  template <class DataType>
  inline void Array3D<DataType>::resize(int const extent_zero,
                                        int const extent_one,
                                        int const extent_two,
                                        DataType const new_value) {
    extent_zero_ = extent_zero;
    extent_one_ = extent_one;
    extent_two_ = extent_two;
    std::size_t const n = extent_zero_ * extent_one_ * extent_two_;
    this->m_.resize(n, new_value);
  }

  template <class DataType>
  inline void Array3D<DataType>::resize(int const extent_zero,
                                        int const extent_one,
                                        int const extent_two,
                                        DataType const *new_array) {
    extent_zero_ = extent_zero;
    extent_one_ = extent_one;
    extent_two_ = extent_two;
    std::size_t const n = extent_zero_ * extent_one_ * extent_two_;
    this->m_.resize(n);
    std::copy(new_array, new_array + n, this->m_.data());
  }

  template <class DataType>
  inline const DataType Array3D<DataType>::operator()(int i, int j, int k) const {
    std::size_t const n = (i * extent_one_ + j) * extent_two_ + k;
    return this->m_[n];
  }

  template <class DataType>
  inline DataType &Array3D<DataType>::operator()(int i, int j, int k) {
    std::size_t const n = (i * extent_one_ + j) * extent_two_ + k;
    return this->m_[n];
  }

  template <class DataType>
  inline DataType const Array3D<DataType>::at(int i, int j, int k) const {
    this->_range_check(i, extent_zero_);
    this->_range_check(j, extent_one_);
    this->_range_check(k, extent_two_);
    std::size_t const n = (i * extent_one_ + j) * extent_two_ + k;
    return this->m_[n];
  }

  template <class DataType>
  inline DataType &Array3D<DataType>::at(int i, int j, int k) {
    this->_range_check(i, extent_zero_);
    this->_range_check(j, extent_one_);
    this->_range_check(k, extent_two_);
    std::size_t const n = (i * extent_one_ + j) * extent_two_ + k;
    return this->m_[n];
  }

  template <class DataType>
  Array3D<DataType>::j_operator::j_operator(Array3D<DataType> &array, int i) : j_array_(array), i_(i) {}

  template <class DataType>
  Array3D<DataType>::j_operator::k_operator::k_operator(Array3D<DataType> &array, int i, int j)
      : k_array_(array), i_(i), j_(j) {}

  template <class DataType>
  inline DataType const Array3D<DataType>::j_operator::k_operator::operator[](int k) const {
    std::size_t const n = (i_ * k_array_.extent_one_ + j_) * k_array_.extent_two_ + k;
    return k_array_.m_[n];
  }

  template <class DataType>
  inline DataType &Array3D<DataType>::j_operator::k_operator::operator[](int k) {
    std::size_t const n = (i_ * k_array_.extent_one_ + j_) * k_array_.extent_two_ + k;
    return k_array_.m_[n];
  }

  template <class DataType>
  inline const typename Array3D<DataType>::j_operator::k_operator Array3D<DataType>::j_operator::operator[](
      int j) const {
    return k_operator(j_array_, i_, j);
  }

  template <class DataType>
  inline typename Array3D<DataType>::j_operator::k_operator Array3D<DataType>::j_operator::operator[](int j) {
    return k_operator(j_array_, i_, j);
  }

  template <class DataType>
  inline const typename Array3D<DataType>::j_operator Array3D<DataType>::operator[](int i) const {
    return j_operator(*this, i);
  }

  template <class DataType>
  inline typename Array3D<DataType>::j_operator Array3D<DataType>::operator[](int i) {
    return j_operator(*this, i);
  }

}  // namespace openKIM

#endif  // HELPER_HPP
