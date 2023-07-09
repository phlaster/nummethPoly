#include "functions.hpp"

Str Buffer::doubleToString(double value, int digits) {
    ostringstream stream;
    stream << scientific << setprecision(digits) << value;
    return stream.str();
}

void Buffer::append(const Str& str) {
    this->data += str;
}

void Buffer::append(int value) {
    this->append(to_string(value)+",");
}

void Buffer::append(double value, bool last) {
    if (!last)
        this->append(doubleToString(value)+",");
    else
        this->append(doubleToString(value)+"\n");
}