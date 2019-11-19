#include "Utils.h"

class Cut {

    public:
        // Constructors
        Cut();
        Cut(std::string name, std::string string_cut);

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

Cut::Cut(std::string name, std::string string_cut) {
    m_name = name;
    m_string_cut = string_cut;
}
