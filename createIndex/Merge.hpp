#ifndef MERGE_HPP
#define MERGE_HPP

#include <utility>

#include <TextPosition.hpp>

/**
 * Type to represent a merge between two bases in a reference structure.
 */
typedef std::pair<std::pair<TextPosition, TextPosition>,bool> Merge;

#endif
