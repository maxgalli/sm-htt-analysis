#ifndef CUTSTRING_H
#define CUTSTRING_H

#include "Utils.h"

class Cut {
    public:
        // Constructors
        Cut();
        Cut(std::string string_cut, std::string name);

        // Getters
        std::string getName() {return m_name;}
        std::string getStringcut() {return m_string_cut;}

        // Setters
        void setName (std::string name) {m_name = name;}
        void setCutstring (std::string& string_cut) {m_string_cut = string_cut;}

    private:
        std::string m_name;
        std::string m_string_cut;
};

#endif
