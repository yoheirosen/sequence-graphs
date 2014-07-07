#include "Mapping.hpp"

Mapping::Mapping(): location(0, 0), is_mapped(false) {
}

Mapping::Mapping(TextPosition location, bool is_mapped, bool is_multimapped): location(location), 
  is_mapped(is_mapped), is_multimapped(is_multimapped) {
}

bool Mapping::operator==(const Mapping& other) const {
  return location == other.location && is_mapped == other.is_mapped && is_multimapped == other.is_multimapped;
}

std::ostream& operator<< (std::ostream& o, Mapping const& mapping) {
    if(mapping.is_mapped) {
        o << "Text " << mapping.location.getText() << " offset " <<
        mapping.location.getOffset();
    } else {
        o << "-----------------";
    }
    return o;
}
