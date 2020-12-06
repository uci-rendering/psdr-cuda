#pragma once

// Object Base Class

class Object {
public:
    virtual ~Object() {}

    virtual std::string type_name() const = 0;

    void log(const char *msg) const {
        std::cout << "[" << type_name() << "] " << msg << std::endl;
    }

    virtual std::string to_string() const {
        std::stringstream oss;
        oss << type_name();
        if ( m_id != "" ) oss << "[id=" << m_id << "]";
        return oss.str();
    }

    std::string m_id = "";
};
